# H&E images are aligned using manually selected landmark pairs from LandmarkPicker.py
# affine transform is performed and then applied to 10x Visium spot coordinates

# load required libraries
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from STalign import STalign

# make plots bigger
plt.rcParams["figure.figsize"] = (12, 10)
# make default project directory
default_dir = Path.cwd()

# expected tissue_positions.csv columns
expected_cols = ["barcode", "in_tissue", "array_row", "array_col",
                 "pxl_row_in_fullres", "pxl_col_in_fullres"]
# Reader of Visium spot data
## making required columns numeric while keeping barcodes as strings
def read_spots(pos_file):
    pos_file = Path(pos_file)
    df = pd.read_csv(pos_file)
    df["barcode"] = df["barcode"].astype(str)
    numeric_cols = ["in_tissue", "array_row", "array_col",
                    "pxl_row_in_fullres","pxl_col_in_fullres"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

# Reader of scale factors
## important is tissue_hires_scalef because has to use high due to landmark choosing
def read_scalefactors(scale_file):
    scale_file = Path(scale_file)
    with open(scale_file, "r") as f:
        scalefactors = json.load(f)
    return scalefactors

# adding x and y coordinates in hires coordinate system (converts full to high resolution)
def add_xy_coordinates(df, scalefactors):
    df = df.copy()
    factor = float(scalefactors["tissue_hires_scalef"])
    df["x"] = df["pxl_col_in_fullres"].astype(float) * factor
    df["y"] = df["pxl_row_in_fullres"].astype(float) * factor
    return df

# Filter of spots to only those in tissue --> needed?
def filter_spots(df):
    df = df.copy()
    df = df[df["in_tissue"].astype(int) == 1].copy()
    df = df.reset_index(drop=True)
    return df

# Spot validation and potential error generation
## checks whether there are in-tissue spots, no duplicate barcodes, only finite x/y coordinates
def validate_spots(df, sample_name):
    if df.empty:
        raise ValueError(f"{sample_name}: no in-tissue spots found")
    if df["barcode"].duplicated().any():
        raise ValueError(f"{sample_name}: duplicate barcodes found")
    if not np.isfinite(df[["x", "y"]].to_numpy()).all():
        raise ValueError(f"{sample_name}: non-finite x/y coordinates found")
    print(f"{sample_name}: {df.shape[0]} in-tissue spots")

## Affine transformation - pre-processing of the data in order to reduce computational burden on STalign
# Reader of manual landmark point CSVs saved from LandmarkPicker.py
def read_landmark_points(points_file):
    points_file = Path(points_file)
    df = pd.read_csv(points_file)
    if not {"y", "x"}.issubset(df.columns):
        raise ValueError(f"{points_file} must contain columns named 'y' and 'x'.")
    points = df[["y", "x"]].astype(float).to_numpy()
    if not np.isfinite(points).all():
        raise ValueError(f"{points_file} contains non-finite landmark coordinates.")
    return points

# Fitter of affine transform from manual landmark pairs (translation, rotation, scaling, shearing)
def fit_affine(points1_yx, points2_yx):
    if points1_yx.shape != points2_yx.shape:
        raise ValueError("Source and reference landmark arrays must have the same shape.")
    if points1_yx.shape[0] < 3:
        raise ValueError("Need at least 3 landmark pairs to fit an affine transform.")
    # converting y,x to x,y for affine fitting
    source_xy = np.column_stack([points1_yx[:, 1], points1_yx[:, 0]])
    reference_xy = np.column_stack([points2_yx[:, 1], points2_yx[:, 0]])
    # adding column of ones for affine translation
    design = np.column_stack([source_xy, np.ones(source_xy.shape[0])])
    # fitting x and y coordinates separately
    coef_x, *_ = np.linalg.lstsq(design, reference_xy[:, 0], rcond=None)
    coef_y, *_ = np.linalg.lstsq(design, reference_xy[:, 1], rcond=None)
    affine = np.vstack([coef_x, coef_y])
    # checking how well landmarks fit after transform
    predicted_xy = design @ affine.T
    residuals = np.linalg.norm(predicted_xy - reference_xy, axis=1)
    return affine, residuals

# applying affine transform to x and y coordinates (transforms source spots)
def apply_affine(x, y, affine):
    xy_hom = np.column_stack([x, y, np.ones(len(x))])
    transformed_xy = xy_hom @ affine.T
    return transformed_xy[:, 0], transformed_xy[:, 1]

# creating output directory and output paths to files
def make_output_paths(args):
    project_dir = Path(args.project_dir)
    if args.outdir is None:
        outdir = (project_dir/ "STalign_outputs"/ f"{args.sample_aligned}_aligned_to_{args.sample_reference}")
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
    return {
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
       "transform_file": outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_{args.alignment_method}_transform.npz"
        )
    }

# printing affine QCs
def print_affine_summary(affine, residuals):
    print(affine)
    print(
        f"Landmark residuals: mean={residuals.mean():.2f}, "
        f"median={np.median(residuals):.2f}, max={residuals.max():.2f}"
    )
    linear_part = affine[:, :2]
    singular_values = np.linalg.svd(linear_part, compute_uv=False)
    determinant = np.linalg.det(linear_part)
    print(f"Affine scale singular values: {singular_values}")
    print(f"Affine determinant: {determinant:.4f}")

## Plots (QC)
# making before-alignment plot
def make_overlay_plot(src, tgt, outpath, title):
    fig, ax = plt.subplots()
    ax.scatter(src["x"], src["y"], s=8, alpha=0.45, label="source")
    ax.scatter(tgt["x"], tgt["y"], s=8, alpha=0.45, label="reference")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# making after-alignment plot
def make_aligned_overlay_plot(src_out, tgt, outpath, title):
    fig, ax = plt.subplots()
    ax.scatter(tgt["x"], tgt["y"], s=8, alpha=0.45, label="reference spots")
    ax.scatter(
        src_out["aligned_x"],
        src_out["aligned_y"],
        s=8,
        alpha=0.45,
        label="source spots aligned",
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# making landmark fit plot (reference, source landmarks after transformation and residual lines between them)
def make_landmark_fit_plot(points1_yx, points2_yx, affine, outpath):
    source_x = points1_yx[:, 1]
    source_y = points1_yx[:, 0]
    reference_x = points2_yx[:, 1]
    reference_y = points2_yx[:, 0]
    aligned_x, aligned_y = apply_affine(source_x, source_y, affine)
    fig, ax = plt.subplots()
    ax.scatter(reference_x, reference_y, s=45, label="reference landmarks")
    ax.scatter(aligned_x, aligned_y, s=45, label="source landmarks after affine")
    # drawing residual lines between aligned source landmarks and reference landmarks
    for x1, y1, x2, y2 in zip(aligned_x, aligned_y, reference_x, reference_y):
        ax.plot([x1, x2], [y1, y2], linewidth=0.8, alpha=0.7)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title("Manual landmark affine fit")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# applying alignment to source spot table and save
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
    # x and y become aligned coordinates in the reference coordinate system, overwrite
    src_out["x"] = src_out["aligned_x"]
    src_out["y"] = src_out["aligned_y"]
    preferred_cols = ["barcode", "x", "y", "original_x", "original_y", "aligned_x", "aligned_y", "in_tissue", "array_row",
                      "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    remaining_cols = [c for c in src_out.columns if c not in preferred_cols]
    src_out = src_out[preferred_cols + remaining_cols]
    return src_out

# saving affine transform and landmark information
def save_transform(transform_file, affine, residuals, points1, points2):
    np.savez_compressed(
        transform_file,
        affine=affine,
        landmark_residuals=residuals,
        points1=points1,
        points2=points2,
        coordinate_system="hires",
    )

##STalign - the module itself building on landmark affine transformation
# reading H&E image for STalign
def read_he_image(image_file):
    image_file = Path(image_file)
    img = plt.imread(image_file)
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1)
    if img.shape[-1] == 4:
        img = img[:, :, :3]
    img = img.astype(float)
    if img.max() > 1:
        img = img / 255.0
    return img

## Preparing H&E image for STalign function
def prepare_stalign_image(img):
    # STalign expects: channels, rows, columns
    img = img.transpose(2, 0, 1)
    img = STalign.normalize(img)
    y = np.arange(img.shape[1]) * 1.0
    x = np.arange(img.shape[2]) * 1.0
    return [y, x], img

# converting torch tensors to numpy
def to_numpy(x):
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)

# running STalign rasterisation of H&E to H&E img
def run_stalign_registration(args, points1, points2):
    if args.image1 is None or args.image2 is None:
        raise ValueError("--image1 and --image2 are required for --alignment_method stalign.")
    device = "cuda:0" if torch.cuda.is_available() else "cpu"

    # reading source and reference H&E images
    img1 = read_he_image(args.image1)
    img2 = read_he_image(args.image2)

    # use function for preparing images for STalign
    xI, I = prepare_stalign_image(img1)
    xJ, J = prepare_stalign_image(img2)

    # no need to flip as LandmarkPicker.py saves y,x, which matches STalign row,column order 
    pointsI = points1
    pointsJ = points2

    # initial affine from landmarks
    L, T = STalign.L_T_from_points(pointsI, pointsJ)

    # running STalign LDDMM registration
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
    return out

# applying STalign transform to source spot table
def align_source_spots_stalign(src, stalign_out):
    # STalign uses y,x (row, column) as mentioned previously
    source_points_yx = np.stack(
        [
            src["y"].to_numpy(),
            src["x"].to_numpy(),
        ],
        axis=1,
    )
    transformed_yx = STalign.transform_points_source_to_target(stalign_out["xv"], stalign_out["v"], stalign_out["A"], source_points_yx)
    transformed_yx = to_numpy(transformed_yx)

    src_out = src.copy()
    src_out["original_x"] = src_out["x"]
    src_out["original_y"] = src_out["y"]

    src_out["aligned_y"] = transformed_yx[:, 0]
    src_out["aligned_x"] = transformed_yx[:, 1]

    # x and y become aligned coordinates in the reference coordinate system
    src_out["x"] = src_out["aligned_x"]
    src_out["y"] = src_out["aligned_y"]

    preferred_cols = ["barcode", "x", "y", "original_x", "original_y", "aligned_x", "aligned_y", "in_tissue", "array_row",
                        "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    remaining_cols = [c for c in src_out.columns if c not in preferred_cols]
    src_out = src_out[preferred_cols + remaining_cols]
    return src_out

# making STalign landmark fit plot
def make_stalign_landmark_fit_plot(points1_yx, points2_yx, stalign_out, outpath):
    transformed_yx = STalign.transform_points_source_to_target(
        stalign_out["xv"],
        stalign_out["v"],
        stalign_out["A"],
        points1_yx,
    )

    transformed_yx = to_numpy(transformed_yx)

    reference_x = points2_yx[:, 1]
    reference_y = points2_yx[:, 0]

    aligned_x = transformed_yx[:, 1]
    aligned_y = transformed_yx[:, 0]

    residuals = np.linalg.norm(transformed_yx - points2_yx, axis=1)

    fig, ax = plt.subplots()
    ax.scatter(reference_x, reference_y, s=45, label="reference landmarks")
    ax.scatter(aligned_x, aligned_y, s=45, label="source landmarks after STalign")

    for x1, y1, x2, y2 in zip(aligned_x, aligned_y, reference_x, reference_y):
        ax.plot([x1, x2], [y1, y2], linewidth=0.8, alpha=0.7)

    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title(
        f"STalign landmark fit / mean residual = {residuals.mean():.2f}px"
    )

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    print(
        f"STalign landmark residuals: mean={residuals.mean():.2f}, "
        f"median={np.median(residuals):.2f}, max={residuals.max():.2f}"
    )

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

# argument parsing
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pos1", required=True, help="Source tissue_positions.csv")
    parser.add_argument("--pos2", required=True, help="Reference tissue_positions.csv")
    parser.add_argument("--scale1", required=True, help="Source scalefactors_json.json")
    parser.add_argument("--scale2", required=True, help="Reference scalefactors_json.json")
    parser.add_argument("--points1", required=True, help="Source manual landmarks CSV with columns y,x")
    parser.add_argument("--points2", required=True, help="Reference manual landmarks CSV with columns y,x")
    parser.add_argument("--project_dir", default=default_dir, type=Path, help="Main project directory where output folders will be created.")
    parser.add_argument("--outdir", default=None, help="Optional custom output directory.")
    parser.add_argument("--outname", default=None)
    parser.add_argument("--sample_aligned", required=True)
    parser.add_argument("--sample_reference", required=True)
    parser.add_argument("--image1", default=None, help="Source tissue_hires_image.png")
    parser.add_argument("--image2", default=None, help="Reference tissue_hires_image.png")
    parser.add_argument("--alignment_method", default="stalign", choices=["affine", "stalign"], help="Alignment method: affine or stalign.")
    parser.add_argument("--niter", default=300, type=int, help="STalign iterations.")
    parser.add_argument("--diffeo_start", default=100, type=int, help="Iteration when nonlinear deformation starts.")
    parser.add_argument("--sigmaM", default=1.5, type=float)
    parser.add_argument("--sigmaB", default=1.0, type=float)
    parser.add_argument("--sigmaA", default=1.1, type=float)
    parser.add_argument("--sigmaP", default=20.0, type=float)
    parser.add_argument("--epV", default=100.0, type=float)
    return parser.parse_args()

# main body using defined functions
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
        stalign_out = run_stalign_registration(args, points1, points2)

        # applying STalign transform to source spots
        src_out = align_source_spots_stalign(src, stalign_out)

        # landmark fit plot
        make_stalign_landmark_fit_plot(points1, points2, stalign_out, paths["landmark_plot"])

        # saving STalign transform
        save_stalign_transform(paths["transform_file"], stalign_out, points1, points2)

    else:
        raise ValueError(f"Other unknown alignment method: {args.alignment_method}")

    # saving aligned coordinate file
    src_out.to_csv(paths["output_csv"], index=False)
    
    # after-alignment plot
    make_aligned_overlay_plot(
        src_out,
        tgt,
        paths["after_plot"],
        title=f"After alignment: {args.sample_aligned} aligned to {args.sample_reference}",
    )

if __name__ == "__main__":
    main()