#// H&E images are aligned using manually selected landmark pairs from LandmarkPicker.py using affine transformation and nonlinear STalign/LDDMM transformation
# there is the possibility to use only affine transformation, but using LDDMM is recommended for better alignment
# the parameters for LDDMM can be adjusted both in the command line and as defaults here in the code
# it is recommended to test varying parameters and check the QC plots to see which parameters work best for your data

# load required libraries
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from STalign import STalign # because the name of the module is the same as of the package 

from qc_plots import (
    to_numpy,
    make_overlay_plot,
    make_aligned_overlay_plot,
    make_stalign_landmark_fit_plot,
    make_stalign_lddmm_diagnostic_plot,
    make_stalign_input_qc_plot,
    make_stalign_initial_affine_qc_plot,
    make_stalign_deformed_image_qc_plot,
    make_stalign_deformation_grid_qc_plot,
    make_stalign_spots_on_target_qc_plot,
    make_stalign_displacement_histogram,
    make_stalign_wm_qc_plots,
    make_qc_summary_panel,
)

# make plots bigger uniformly 
plt.rcParams["figure.figsize"] = (12, 10)

# make default project directory
default_dir = Path.cwd()

#// function reading Visium spot data
# making required columns numeric while keeping barcodes as strings
def read_spots(pos_file):
    pos_file = Path(pos_file)
    df = pd.read_csv(pos_file)
    df["barcode"] = df["barcode"].astype(str)
    numeric_cols = ["in_tissue", "array_row", "array_col",
                    "pxl_row_in_fullres","pxl_col_in_fullres"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="raise")
    return df

#// function reading scale factors
# important is tissue_hires_scalef because has to use high due to landmark choosing in LandmarkPicker.py
def read_scalefactors(scale_file):
    scale_file = Path(scale_file)
    with open(scale_file, "r") as f:
        scalefactors = json.load(f)
    return scalefactors

#// function adding x and y coordinates in hires coordinate system
# converts full to high resolution for less computational burden
def add_xy_coordinates(df, scalefactors):
    factor = float(scalefactors["tissue_hires_scalef"])
    df["x"] = df["pxl_col_in_fullres"].astype(float) * factor
    df["y"] = df["pxl_row_in_fullres"].astype(float) * factor
    return df

#// function filtering spots to only those in tissue 
def filter_spots(df):
    df = df[df["in_tissue"].astype(int) == 1].copy()
    df = df.reset_index(drop=True)
    return df

#// function validating spots
## checks whether there are in-tissue spots, no duplicate barcodes, only finite x/y coordinates
def validate_spots(df, sample_name):
    if df.empty:
        raise ValueError(f"{sample_name}: no in-tissue spots found")
    if df["barcode"].duplicated().any():
        raise ValueError(f"{sample_name}: duplicate barcodes found")
    if not np.isfinite(df[["x", "y"]].to_numpy()).all():
        raise ValueError(f"{sample_name}: non-finite x/y coordinates found")
    print(f"{sample_name}: {df.shape[0]} in-tissue spots") #  is for rows who are barcodes and thus spots 

## Affine transformation
#// function reading manual landmark points from CSVs saved from LandmarkPicker.py
# returns numpy array with y,x coordinates (this order is required for further analysis)
def read_landmark_points(points_file):
    points_file = Path(points_file)
    df = pd.read_csv(points_file)
    if not {"y", "x"}.issubset(df.columns):
        raise ValueError(f"{points_file} must contain columns named 'y' and 'x'.")
    points = df[["y", "x"]].astype(float).to_numpy()
    if not np.isfinite(points).all():
        raise ValueError(f"{points_file} contains non-finite landmark coordinates.")
    return points

#// function fitting affine transformation from source to reference landmarks
# returns affine matrix and residuals of the fit
# first converts y, x to x, y for fitting, then adds column of ones to allow transformation
# fits x and y coordinates separately using least squares
# returns matrix and residuals for the QC of the fit 
def fit_affine(points1_yx, points2_yx):
    if points1_yx.shape != points2_yx.shape:
        raise ValueError("Source and reference landmark arrays must have the same shape.")

    source_xy = np.column_stack([points1_yx[:, 1], points1_yx[:, 0]])
    reference_xy = np.column_stack([points2_yx[:, 1], points2_yx[:, 0]])
  
    design = np.column_stack([source_xy, np.ones(source_xy.shape[0])])
    coef_x, *_ = np.linalg.lstsq(design, reference_xy[:, 0], rcond=None)
    coef_y, *_ = np.linalg.lstsq(design, reference_xy[:, 1], rcond=None)
    affine = np.vstack([coef_x, coef_y])
  
    predicted_xy = design @ affine.T # applied to source landmarks coords and compares
    residuals = np.linalg.norm(predicted_xy - reference_xy, axis=1)
    return affine, residuals

#// function applying affine transformation to x and y coordinates (landmarks and source spots)
def apply_affine(x, y, affine):
    xy_hom = np.column_stack([x, y, np.ones(len(x))])
    transformed_xy = xy_hom @ affine.T 
    return transformed_xy[:, 0], transformed_xy[:, 1]

#// function creating output directory and output paths to files
def make_output_paths(args):
    project_dir = Path(args.project_dir)
    if args.outdir is None:
        outdir = (project_dir/ "STalign_outputs"/ f"{args.sample_aligned}_aligned_to_{args.sample_reference}") # could eb changed to be cwd
    else:
        outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.outname is None:
        outname = (
            f"{args.sample_aligned}_aligned_to_"
            f"{args.sample_reference}_barcodes.csv"
        )
    else:
        outname = args.outname 
    d = {
        "outdir": outdir,
        "output_csv": outdir / outname,
        "before_plot": outdir / (
            f"{args.sample_aligned}_vs_{args.sample_reference}_before_alignment.png"
        ),
        "after_plot": outdir / (
            f"{args.sample_aligned}_aligned_to_{args.sample_reference}_after_alignment.png"
        ),
        "landmark_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_manual_landmark_fit.png"
        ),
        "stalign_input_qc_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_inputs_landmarks.png"
        ),
        "stalign_initial_affine_qc_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_initial_affine.png"
        ),
        "stalign_deformed_image_qc_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_deformed_source_vs_target.png"
        ),
        "stalign_deformation_grid_qc_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_deformation_grid.png"
        ),
        "stalign_spots_on_target_qc_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_aligned_spots_on_target.png"
        ),
        "stalign_displacement_hist_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_spot_displacement_histogram.png"
        ),
        "stalign_wm_spot_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_WM_values_on_spots.png"
        ),
        "stalign_wm_hist_plot": outdir / (
            f"{args.sample_aligned}_to_{args.sample_reference}_WM_values_histogram.png"
        ),
       "transform_file": outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_{args.alignment_method}_transform.npz"
        ),
        "stalign_lddmm_diagnostic_plot": outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_LDDMM_diagnostic.png"
        ),
        "qc_summary_plot": outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_QC_summary_panel.png"
        ),
    }
    return d

#// function printing affine matrix and residuals of the fit
def print_affine_summary(affine, residuals):
    print(affine)
    print(
        f"Landmark residuals: mean={residuals.mean():.2f}, "
        f"median={np.median(residuals):.2f}, max={residuals.max():.2f}"
    )
    linear_part = affine[:, :2]
    singular_values = np.linalg.svd(linear_part, compute_uv=False) # stretching, should remain similar on both sides 
    determinant = np.linalg.det(linear_part) # should be positive, if negative then the image is flipped
    print(f"Affine scale singular values: {singular_values}")
    print(f"Affine determinant: {determinant:.2f}")

#// function applying alignment to source spot table 
def align_source_spots(src, affine):
    aligned_x, aligned_y = apply_affine(
        src["x"].to_numpy(),
        src["y"].to_numpy(),
        affine,
    )
    src_out = src.copy()
    src_out["original_x"] = src_out["x"]
    src_out["original_y"] = src_out["y"]
    src_out["aligned_x"] = aligned_x
    src_out["aligned_y"] = aligned_y

    src_out["x"] = src_out["aligned_x"]
    src_out["y"] = src_out["aligned_y"]
    preferred_cols = ["barcode", "x", "y", "original_x", "original_y", "aligned_x", "aligned_y", "in_tissue", "array_row",
                      "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    src_out = src_out[preferred_cols]
    return src_out

#// function making landmark fit plot (reference, source landmarks after transformation and residual lines between them)
# draws residual lines between aligned source landmarks and reference landmarks
def make_landmark_fit_plot(points1_yx, points2_yx, affine, outpath):
    source_x = points1_yx[:, 1]
    source_y = points1_yx[:, 0]
    reference_x = points2_yx[:, 1]
    reference_y = points2_yx[:, 0]
    aligned_x, aligned_y = apply_affine(source_x, source_y, affine)
    fig, ax = plt.subplots()
    ax.scatter(reference_x, reference_y, s=45, label="reference landmarks")
    ax.scatter(aligned_x, aligned_y, s=45, label="source landmarks after affine")

    for x1, y1, x2, y2 in zip(aligned_x, aligned_y, reference_x, reference_y):
        ax.plot([x1, x2], [y1, y2], linewidth=0.8, alpha=0.7)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title("Manual landmark affine fit")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

#// function saving affine transform and landmark information
def save_transform(transform_file, affine, residuals, points1, points2):
    np.savez_compressed(
        transform_file,
        affine=affine,
        landmark_residuals=residuals,
        points1=points1,
        points2=points2,
        coordinate_system="hires",
    )

##STalign - nonlinear transformation (LDDMM)
# reading H&E image for STalign and checking the number of channels and rescaling if not on 0-1 scale 
def read_he_image(image_file):
    image_file = Path(image_file)
    img = plt.imread(image_file)
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1) # converted greyscale 2D to 3D RGB
    if img.shape[-1] == 4:
        img = img[:, :, :3] # 4th channel is not needed, removed to 3D RGB
    img = img.astype(float)
    if img.max() > 1:
        img = img / 255.0 # might be 8bit and thus require normalisation
    return img

## Preparing H&E images
#// function preparing H&E image for STalign function
# converts channel-last image to channel-first and normalizes to zero mean and unit variance
# generates y and x axes for the image in STalign format (y is rows, x is columns)
# converts to float 
def prepare_stalign_image(img):
    img = img.transpose(2, 0, 1)
    img = STalign.normalize(img)
    y = np.arange(img.shape[1]) * 1.0 
    x = np.arange(img.shape[2]) * 1.0
    return [y, x], img

#// function checking image input for STalign
# checks that the images are present
# chooses computation device (GPU if available, otherwise CPU)
def run_stalign_registration(args, points1, points2):
    if args.image1 is None or args.image2 is None:
        raise ValueError("--image1 and --image2 are required for --alignment_method stalign.")
    if torch.cuda.is_available():
        device = "cuda:0"
    else:
        device = "cpu"

    ## Affine transformation

    # reading source and reference H&E images
    img1 = read_he_image(args.image1)
    img2 = read_he_image(args.image2)

    # xI/J are coordinate arrays for the images (in y and x order)
    # I/J are actual images in STalign format (normalised and channel-first)
    xI, I = prepare_stalign_image(img1)
    xJ, J = prepare_stalign_image(img2)

    # pointsI/J are the landmark points in y,x order (STalign requires this order, already saved as such)
    pointsI = points1
    pointsJ = points2

    # initial affine from landmarks
    # L are landmark points in source image, T are landmark points in reference image (target)
    L, T = STalign.L_T_from_points(pointsI, pointsJ)

    # converts L and T to torch tensors (needed for STalign), then calculates initial affine transformation A
    A_init = STalign.to_A(
        torch.as_tensor(to_numpy(L), dtype=torch.float64, device=device),
        torch.as_tensor(to_numpy(T), dtype=torch.float64, device=device),
    )

    ## Non-linear LDDMM transformation

    # parameters for nonlinear STalign/LDDMM (can be adjusted here or in the command line)
    params = {
        "L": L,
        "T": T,
        "pointsI": pointsI,
        "pointsJ": pointsJ,
        "niter": args.niter,
        "diffeo_start": args.diffeo_start,
        "device": device,
        "sigmaM": args.sigmaM,
        "sigmaB": args.sigmaB,
        "sigmaA": args.sigmaA,
        "sigmaP": args.sigmaP,
        "epV": args.epV,
    }
    out = STalign.LDDMM(xI, I, xJ, J, **params)
    print("STalign output keys:", list(out.keys()))

    # diagnostic information about STalign weight maps (QC)
    # WM - matching tissue weight map, WB - background weight map, WA - artefact weight map
    for key in ["WM", "WB", "WA"]:
     if key in out:
        arr = to_numpy(out[key])
        arr = np.squeeze(arr)
        print(f"{key} shape:", arr.shape)
        print(f"{key} min:", np.nanmin(arr))
        print(f"{key} max:", np.nanmax(arr))
        print(f"{key} mean:", np.nanmean(arr))
        print(f"{key} nonzero fraction:", np.mean(arr > 0))

    stalign_data = {
        "xI": xI,
        "I": I,
        "xJ": xJ,
        "J": J,
        "A_init": A_init,
        "device": device,
    }
    return out, stalign_data

#// function transforming source spots to reference coordinate system using STalign output
def align_source_spots_stalign(src, stalign_out):
    source_points_yx = src[["y", "x"]].to_numpy(dtype=float)
    transformed_yx = STalign.transform_points_source_to_target(stalign_out["xv"], stalign_out["v"], stalign_out["A"], source_points_yx)
    transformed_yx = to_numpy(transformed_yx)

# adding transformed coordinates to the df while preserving original coordinates  
    src_out = src.copy()
    src_out["original_x"] = src_out["x"]
    src_out["original_y"] = src_out["y"]
    src_out["aligned_y"] = transformed_yx[:, 0] 
    src_out["aligned_x"] = transformed_yx[:, 1]
    src_out["x"] = src_out["aligned_x"]
    src_out["y"] = src_out["aligned_y"]

    preferred_cols = ["barcode", "x", "y", "original_x", "original_y", "aligned_x", "aligned_y", "in_tissue", "array_row",
                        "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    
    remaining_cols = [c for c in src_out.columns if c not in preferred_cols]
    src_out = src_out[preferred_cols + remaining_cols]
    return src_out

#// function adding spot-level STalign QC metrics for both affine and LDDMM transformation
def stalign_qc(src_out, stalign_out, stalign_data):
    src_out = src_out.copy()

    # total movement from original source coordinates to final target-space coordinates (both affine and LDDMM)
    src_out["stalign_total_displacement"] = np.sqrt(
        (src_out["aligned_x"] - src_out["original_x"]) ** 2
        + (src_out["aligned_y"] - src_out["original_y"]) ** 2
    )
    src_out["stalign_displacement"] = src_out["stalign_total_displacement"]

    # movement added by nonlinear LDDMM beyond the initial landmark affine (on top of affine, check for high values)
    if "A_init" in stalign_data:
        A_init_np = to_numpy(stalign_data["A_init"])
        source_hom = np.vstack(
            [
                src_out["original_y"].to_numpy(),
                src_out["original_x"].to_numpy(),
                np.ones(src_out.shape[0]),
            ]
        )
        affine_spots = (A_init_np @ source_hom).T 
        src_out["stalign_affine_y"] = affine_spots[:, 0]
        src_out["stalign_affine_x"] = affine_spots[:, 1]
        src_out["stalign_nonlinear_displacement"] = np.sqrt(
            (src_out["aligned_x"] - src_out["stalign_affine_x"]) ** 2
            + (src_out["aligned_y"] - src_out["stalign_affine_y"]) ** 2
        )

    if "WM" in stalign_out:
        try:
            spot_yx = np.stack(
                [
                    src_out["aligned_y"].to_numpy(),
                    src_out["aligned_x"].to_numpy(),
                ],
                axis=1,
            )
            wm = stalign_out["WM"]
            if not isinstance(wm, torch.Tensor):
                wm = torch.as_tensor(wm)
            spot_yx_tensor = torch.as_tensor(
                spot_yx,
                dtype=torch.float64,
                device=wm.device,
            )
            wm_values = STalign.interp(
                stalign_data["xJ"],
                wm[None].float(),
                spot_yx_tensor[None].permute(-1, 0, 1).float(),
            )
            src_out["stalign_WM_value"] = to_numpy(wm_values).squeeze()
        except Exception as e:
            print(f"Warning: could not sample WM values at aligned spot positions: {e}")
    return src_out

# saving STalign transform
def save_stalign_transform(transform_file, stalign_out, points1, points2):
    np.savez_compressed(
        transform_file,
        A=to_numpy(stalign_out["A"]),
        v=to_numpy(stalign_out["v"]),
        xv0=to_numpy(stalign_out["xv"][0]),
        xv1=to_numpy(stalign_out["xv"][1]),
        points1=points1,
        points2=points2,
        coordinate_system="hires",
        method="stalign",
    )

# argument parsing (using argparse)
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pos1", required=True, help="Source tissue_positions.csv")
    parser.add_argument("--pos2", required=True, help="Reference tissue_positions.csv")
    parser.add_argument("--scale1", required=True, help="Source scalefactors_json.json")
    parser.add_argument("--scale2", required=True, help="Reference scalefactors_json.json")
    parser.add_argument("--points1", required=True, help="Source manual landmarks CSV with columns y,x")
    parser.add_argument("--points2", required=True, help="Reference manual landmarks CSV with columns y,x")
    parser.add_argument("--project_dir", default=default_dir, type=Path, help="Main project directory where output folders will be created.")
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--outname", default=None)
    parser.add_argument("--sample_aligned", required=True)
    parser.add_argument("--sample_reference", required=True)
    parser.add_argument("--image1", default=None, help="Source tissue_hires_image.png")
    parser.add_argument("--image2", default=None, help="Reference tissue_hires_image.png")
    parser.add_argument("--alignment_method", default="stalign", choices=["affine", "stalign"], help="Alignment method: affine or stalign.")
    parser.add_argument("--niter", default=500, type=int, help="STalign iterations.")
    parser.add_argument("--diffeo_start", default=100, type=int, help="Iteration when nonlinear deformation starts.")
    parser.add_argument("--sigmaM", default=0.18, type=float)
    parser.add_argument("--sigmaB", default=0.18, type=float)
    parser.add_argument("--sigmaA", default=0.18, type=float)
    parser.add_argument("--sigmaP", default=2e-1, type=float)
    parser.add_argument("--epV", default=5e1, type=float)
    parser.add_argument("--skip_stalign_qc", action="store_true", help="Skip extra STalign QC plots.")
    parser.add_argument("--stalign_grid_levels", default=20, type=int, help="Number of contour levels for the deformation-grid QC plot.")
    return parser.parse_args() 

def main():
    # parsing of arguments and path making
    args = parse_args()
    paths = make_output_paths(args)

    # reading source sample
    src_raw = read_spots(args.pos1)
    src_sf = read_scalefactors(args.scale1)

    # reading reference sample
    tgt_raw = read_spots(args.pos2)
    tgt_sf = read_scalefactors(args.scale2)

    # adding coordinates
    src = add_xy_coordinates(src_raw, src_sf)
    tgt = add_xy_coordinates(tgt_raw, tgt_sf)

    # filtering to in-tissue spots
    src = filter_spots(src)
    tgt = filter_spots(tgt)

    # validation
    validate_spots(src, args.sample_aligned)
    validate_spots(tgt, args.sample_reference)

    # reading manual landmarks
    points1 = read_landmark_points(args.points1)
    points2 = read_landmark_points(args.points2)

    # before-alignment plot
    make_overlay_plot(src, tgt, paths["before_plot"], title=f"Before alignment: {args.sample_aligned} vs {args.sample_reference}")

## Affine version 
    if args.alignment_method == "affine":
        # fitting affine transform
        affine, residuals = fit_affine(points1, points2)
        print_affine_summary(affine, residuals)

        # applying affine transform to source spots
        src_out = align_source_spots(src, affine)

        # landmark fit plot
        make_landmark_fit_plot(points1, points2, affine, paths["landmark_plot"])

        # saving affine transform
        save_transform(paths["transform_file"], affine, residuals, points1, points2)

## STalign version
    elif args.alignment_method == "stalign":
        # running STalign H&E-to-H&E registration
        stalign_out, stalign_data = run_stalign_registration(args, points1, points2)

        # applying STalign transform to source spots
        src_out = align_source_spots_stalign(src, stalign_out)

        # add spot-level QC metrics before saving the output table
        src_out = stalign_qc(
            src_out,
            stalign_out,
            stalign_data,
        )

        # landmark fit plot
        make_stalign_landmark_fit_plot(points1, points2, stalign_out, paths["landmark_plot"])

        # STalign QC plots
        if not args.skip_stalign_qc:
            make_stalign_input_qc_plot(
                stalign_data,
                points1,
                points2,
                paths["stalign_input_qc_plot"],
            )
            make_stalign_initial_affine_qc_plot(
                stalign_data,
                points1,
                points2,
                paths["stalign_initial_affine_qc_plot"],
            )
            make_stalign_deformed_image_qc_plot(
                stalign_data,
                stalign_out,
                paths["stalign_deformed_image_qc_plot"],
            )
            make_stalign_deformation_grid_qc_plot(
                stalign_data,
                stalign_out,
                points1,
                points2,
                paths["stalign_deformation_grid_qc_plot"],
                grid_levels=args.stalign_grid_levels,
            )
            make_stalign_spots_on_target_qc_plot(
                src_out,
                tgt,
                stalign_data,
                paths["stalign_spots_on_target_qc_plot"],
                title=f"STalign spot QC: {args.sample_aligned} aligned to {args.sample_reference}",
            )
            make_stalign_displacement_histogram(
                src_out,
                paths["stalign_displacement_hist_plot"],
            )
            make_stalign_wm_qc_plots(
                src_out,
                paths["stalign_wm_spot_plot"],
                paths["stalign_wm_hist_plot"],
            )
            make_stalign_lddmm_diagnostic_plot(
                stalign_data,
                stalign_out,
                paths["stalign_lddmm_diagnostic_plot"],
            )
        # saving STalign transform
        save_stalign_transform(paths["transform_file"], stalign_out, points1, points2)

    else:
        raise ValueError(f"Other unknown alignment method: {args.alignment_method}")

    # saving aligned coordinate file
    src_out.to_csv(paths["output_csv"], index=False)
    
    # after-alignment plot
    make_aligned_overlay_plot(src_out, tgt, paths["after_plot"], title=f"After alignment: {args.sample_aligned} aligned to {args.sample_reference}")
    # summary panel 
    make_qc_summary_panel(paths, args.alignment_method, paths["qc_summary_plot"], title=f"QC summary: {args.sample_aligned} aligned to {args.sample_reference}")

if __name__ == "__main__":
    main()