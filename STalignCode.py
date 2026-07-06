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
from STalign import STalign # because the name of the module is the same as of the package 

# make plots bigger uniformly 
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
        df[col] = pd.to_numeric(df[col], errors="raise")
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
    factor = float(scalefactors["tissue_hires_scalef"])
    df["x"] = df["pxl_col_in_fullres"].astype(float) * factor
    df["y"] = df["pxl_row_in_fullres"].astype(float) * factor
    return df

# Filter of spots to only those in tissue --> good to have
def filter_spots(df):
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
    print(f"{sample_name}: {df.shape[0]} in-tissue spots") #  is for rows who are barcodes and thus spots 

## Affine transformation - reader of pre-processing of the data in order to reduce computational burden on STalign
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
    predicted_xy = design @ affine.T # applied to source landmarks coords and compares
    residuals = np.linalg.norm(predicted_xy - reference_xy, axis=1)
    return affine, residuals

# applying affine transform to x and y coordinates (transforms source spots)
def apply_affine(x, y, affine):
    xy_hom = np.column_stack([x, y, np.ones(len(x))])
    transformed_xy = xy_hom @ affine.T # applied to coordinate matrix  
    return transformed_xy[:, 0], transformed_xy[:, 1]

# creating output directory and output paths to files
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
    d = { # conventional way
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
        # Extra QC plots for STalign/LDDMM. These are only produced when
        # --alignment_method stalign is used.
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
            f"{args.sample_aligned}_to_{args.sample_reference}__WM_values_histogram.png"
        ),
       "transform_file": outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_{args.alignment_method}_transform.npz"
        ),
        "stalign_lddmm_diagnostic_plot": outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_LDDMM_diagnostic.png"
),
    }
    return d

# printing affine QCs
def print_affine_summary(affine, residuals):
    print(affine)
    print(
        f"Landmark residuals: mean={residuals.mean():.2f}, "
        f"median={np.median(residuals):.2f}, max={residuals.max():.2f}"
    )
    linear_part = affine[:, :2]
    singular_values = np.linalg.svd(linear_part, compute_uv=False) # stretching, should remain similar on both sides, no uv has to be output 
    determinant = np.linalg.det(linear_part) # positive means no flipping (should be positive, not negative and flipped and not 0 which would be line), ideally around 1 to have preserved orientation
    print(f"Affine scale singular values: {singular_values}")
    print(f"Affine determinant: {determinant:.2f}")

## Plots (QC)
# making before-alignment plot
def make_overlay_plot(src, tgt, outpath, title):
    fig, ax = plt.subplots()
    ax.scatter(src["x"], src["y"], s=8, alpha=0.45, label="source")
    ax.scatter(tgt["x"], tgt["y"], s=8, alpha=0.45, label="reference")
    ax.set_aspect("equal")
    ax.invert_yaxis() # key
    ax.legend()
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# making after-alignment plot
def make_aligned_overlay_plot(src_out, tgt, outpath, title):
    fig, ax = plt.subplots()
    ax.scatter(tgt["x"], tgt["y"], s=8, alpha=0.45, label="reference spots")
    ax.scatter(src_out["aligned_x"], src_out["aligned_y"], s=8, alpha=0.45,label="source spots aligned")
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

# applying alignment to source spot table and save ## perhaps move earlier 
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
    src_out = src_out[preferred_cols]
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
# reading H&E image for STalign and checking the number of channels and rescaling if not on 0-1 scale 
def read_he_image(image_file):
    image_file = Path(image_file)
    img = plt.imread(image_file)
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1) # converted greyscale 2D to 3D RGB
    if img.shape[-1] == 4:
        img = img[:, :, :3] # 4th channel with transparency is not needed, only RGB, removed 
    img = img.astype(float)
    if img.max() > 1:
        img = img / 255.0 # potential present cmmon 255 format where 0 i bcgrnd and 25 is white, normalise if needed (no other included --> might)
    return img

## Preparing H&E image for STalign function
def prepare_stalign_image(img):
    # STalign expects: channels, rows, columns (normally channels last)
    img = img.transpose(2, 0, 1)
    img = STalign.normalize(img)
    y = np.arange(img.shape[1]) * 1.0 # generates axes and converts to floats
    x = np.arange(img.shape[2]) * 1.0
    return [y, x], img

# converting torch tensors to numpy
def to_numpy(x):
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)

# helper for plotting STalign images regardless of whether they are torch tensors,
# numpy arrays, channel-first, or channel-last arrays
def image_for_plot(img):
    arr = to_numpy(img)
    if arr.ndim == 3 and arr.shape[0] in (1, 3, 4):
        arr = np.moveaxis(arr, 0, -1)
    if arr.ndim == 3 and arr.shape[-1] == 1:
        arr = arr[:, :, 0]
    arr = arr.astype(float)
    finite = np.isfinite(arr)
    if finite.any():
        lo = np.nanmin(arr[finite])
        hi = np.nanmax(arr[finite])
        if hi > lo:
            arr = (arr - lo) / (hi - lo)
    return arr

# helper for consistent image extents from STalign y,x coordinate axes
def stalign_extent(x_axes):
    try:
        return STalign.extent_from_x(tuple(x_axes))
    except Exception:
        y_axis, x_axis = x_axes
        return [x_axis[0], x_axis[-1], y_axis[-1], y_axis[0]]

# robust contour helper so deformation-grid plotting does not crash when values are flat
def add_contours(ax, x_axis, y_axis, z, n_levels=20, **kwargs):
    z = to_numpy(z)
    finite = np.isfinite(z)
    if not finite.any():
        return
    lo, hi = np.nanpercentile(z[finite], [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return
    levels = np.linspace(lo, hi, n_levels)
    X, Y = np.meshgrid(x_axis, y_axis)
    ax.contour(X, Y, z, levels=levels, linewidths=0.5, alpha=0.7, **kwargs)

# running STalign registration of H&E to H&E img
def run_stalign_registration(args, points1, points2):
    if args.image1 is None or args.image2 is None:
        raise ValueError("--image1 and --image2 are required for --alignment_method stalign.")
    if torch.cuda.is_available():
        device = "cuda:0"
    else:
        device = "cpu"

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

    # save the initial affine matrix as a plotting QC before nonlinear LDDMM
    # The notebook visualizes this step before deciding whether LDDMM is needed.
    A_init = STalign.to_A(
        torch.as_tensor(to_numpy(L), dtype=torch.float64, device=device),
        torch.as_tensor(to_numpy(T), dtype=torch.float64, device=device),
    )

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
    print("STalign output keys:", list(out.keys()))

    for key in ["WM", "WB", "WA"]:
     if key in out:
        arr = to_numpy(out[key])
        arr = np.squeeze(arr)
        print(f"{key} shape:", arr.shape)
        print(f"{key} min:", np.nanmin(arr))
        print(f"{key} max:", np.nanmax(arr))
        print(f"{key} mean:", np.nanmean(arr))
        print(f"{key} nonzero fraction:", np.mean(arr > 0))
    # keep image inputs/axes because the extra QC plots need them later
    stalign_data = {
        "xI": xI,
        "I": I,
        "xJ": xJ,
        "J": J,
        "A_init": A_init,
        "device": device,
    }
    return out, stalign_data

# applying STalign transform to source spot table
def align_source_spots_stalign(src, stalign_out):
    # STalign uses y,x (row, column) as mentioned previously
    source_points_yx = np.stack( # might be unnecessary 
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

    src_out["aligned_y"] = transformed_yx[:, 0] # 0 is y remember 
    src_out["aligned_x"] = transformed_yx[:, 1]

    # x and y become aligned coordinates in the reference coordinate system
    src_out["x"] = src_out["aligned_x"]
    src_out["y"] = src_out["aligned_y"]

    preferred_cols = ["barcode", "x", "y", "original_x", "original_y", "aligned_x", "aligned_y", "in_tissue", "array_row",
                        "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    remaining_cols = [c for c in src_out.columns if c not in preferred_cols]
    src_out = src_out[preferred_cols + remaining_cols]
    return src_out

# making STalign fit plot
def make_stalign_landmark_fit_plot(points1_yx, points2_yx, stalign_out, outpath):
    transformed_yx = STalign.transform_points_source_to_target(
        stalign_out["xv"],
        stalign_out["v"],
        stalign_out["A"],
        points1_yx,
    )
## where is this?
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
####
def make_stalign_lddmm_diagnostic_plot(stalign_data, stalign_out, outpath):
    xI = stalign_data["xI"]
    I = stalign_data["I"]
    xJ = stalign_data["xJ"]
    J = stalign_data["J"]

    xv = stalign_out["xv"]
    v = stalign_out["v"]
    A = stalign_out["A"]
    WM = stalign_out.get("WM", None)



    # final transformed source image
    phiI = STalign.transform_image_source_to_target(
        xv,
        v,
        A,
        xI,
        I,
        xJ,
    )

    # deformation transform on target grid
    phii = STalign.build_transform(
        xv,
        v,
        A,
        XJ=xJ,
        direction="b",
    )

    I_plot = image_for_plot(I)
    J_plot = image_for_plot(J)
    phiI_plot = image_for_plot(phiI)

    yJ, xJ_axis = xJ
    extentJ = stalign_extent(xJ)

    # crude residual image: target - transformed source
    if phiI_plot.ndim == 3 and J_plot.ndim == 3:
        error_img = np.mean(np.abs(J_plot - phiI_plot), axis=2)
    else:
        error_img = np.abs(np.asarray(J_plot) - np.asarray(phiI_plot))

    fig, ax = plt.subplots(2, 3, figsize=(18, 11))

    ax[0, 0].imshow(phiI_plot, extent=extentJ)
    ax[0, 0].set_title("Space/contrast transformed source")

    ax[0, 1].imshow(J_plot, extent=extentJ)
    ax[0, 1].set_title("Target")

    ax[0, 2].imshow(error_img, extent=extentJ)
    ax[0, 2].set_title("Residual error")

    phii_np = to_numpy(phii)
    if phii_np.ndim >= 3:
        # approximate deformation magnitude
        yy, xx = np.meshgrid(yJ, xJ_axis, indexing="ij")
        deformation_mag = np.sqrt(
            (phii_np[..., 0] - yy) ** 2 +
            (phii_np[..., 1] - xx) ** 2
        )
        im = ax[1, 0].imshow(deformation_mag, extent=extentJ)
        ax[1, 0].set_title("Deformation magnitude")
        fig.colorbar(im, ax=ax[1, 0], fraction=0.046, pad=0.04)

    if WM is not None:
        WM_np = to_numpy(WM)
        im = ax[1, 1].imshow(WM_np, extent=extentJ)
        ax[1, 1].set_title("STalign weights / WM")
        fig.colorbar(im, ax=ax[1, 1], fraction=0.046, pad=0.04)
    else:
        ax[1, 1].axis("off")
        ax[1, 1].set_title("No WM found")

    ax[1, 2].imshow(phiI_plot, extent=extentJ)
    if phii_np.ndim >= 3:
        add_contours(ax[1, 2], xJ_axis, yJ, phii_np[..., 0], n_levels=20)
        add_contours(ax[1, 2], xJ_axis, yJ, phii_np[..., 1], n_levels=20)
    ax[1, 2].set_title("Transformed source with deformation grid")

    for a in ax.ravel():
        a.set_aspect("equal")
        a.invert_yaxis()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

    #### 
# STalign QC metrics to check
# For Visium, each row is a spot (not a sc) so WM values are attached to spots as QC information, but the spots are not filtered by default
def stalign_qc(src_out, stalign_out, stalign_data, wm_threshold=None):
    src_out = src_out.copy()

    # Total movement from original source coordinates to final target-space coordinates (both affine and LDDMM)
    src_out["stalign_total_displacement"] = np.sqrt(
        (src_out["aligned_x"] - src_out["original_x"]) ** 2
        + (src_out["aligned_y"] - src_out["original_y"]) ** 2
    )
    src_out["stalign_displacement"] = src_out["stalign_total_displacement"]

    # Movement added by nonlinear LDDMM beyond the initial landmark affine (on top of affine, check for high values)
    if "A_init" in stalign_data:
        A_init_np = to_numpy(stalign_data["A_init"])
        source_hom = np.vstack(
            [
                src_out["original_y"].to_numpy(),
                src_out["original_x"].to_numpy(),
                np.ones(src_out.shape[0]),
            ]
        )
        affine_spots = (A_init_np @ source_hom).T # T?
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
            src_out["stalign_WM_value"] = to_numpy(wm_values[0, 0])
            if wm_threshold is not None:
                src_out["stalign_WM_pass"] = src_out["stalign_WM_value"] >= wm_threshold
        except Exception as e:
            print(f"Could not compute STalign WM values for spots: {e}")
    return src_out

# QC 1: visualize the exact H&E images (ovrelapping) and landmarks given to STalign
def make_stalign_input_qc_plot(stalign_data, points1_yx, points2_yx, outpath):
    xI = stalign_data["xI"]
    I = stalign_data["I"]
    xJ = stalign_data["xJ"]
    J = stalign_data["J"]

    fig, ax = plt.subplots(1, 2, figsize=(16, 8))
    ax[0].imshow(image_for_plot(I), extent=stalign_extent(xI))
    ax[0].scatter(points1_yx[:, 1], points1_yx[:, 0], s=35, label="source landmarks")
    ax[0].set_title("STalign source H&E with landmarks")
    ax[0].set_aspect("equal")
    ax[0].invert_yaxis()
    ax[0].legend()

    ax[1].imshow(image_for_plot(J), extent=stalign_extent(xJ))
    ax[1].scatter(points2_yx[:, 1], points2_yx[:, 0], s=35, label="reference landmarks")
    ax[1].set_title("STalign reference H&E with landmarks")
    ax[1].set_aspect("equal")
    ax[1].invert_yaxis()
    ax[1].legend()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# QC 2: affine alignment before LDDMM
def make_stalign_initial_affine_qc_plot(stalign_data, points1_yx, points2_yx, outpath):
    xI = stalign_data["xI"]
    I = stalign_data["I"]
    xJ = stalign_data["xJ"]
    J = stalign_data["J"]
    A_init = stalign_data["A_init"]

    affine_I = STalign.transform_image_source_with_A(A_init, xI, I, xJ)
    # applying only the affine transformation which accepts the coordinates in the matrix as y and x, therefore keep the alignment 
    A_init_np = to_numpy(A_init)
    source_hom = np.vstack(
        [
            points1_yx[:, 0],
            points1_yx[:, 1],
            np.ones(points1_yx.shape[0]),
        ]
    )
    transformed_landmarks = (A_init_np @ source_hom).T

    fig, ax = plt.subplots(1, 2, figsize=(16, 8))
    ax[0].imshow(image_for_plot(affine_I), extent=stalign_extent(xJ))
    ax[0].scatter(transformed_landmarks[:, 1], transformed_landmarks[:, 0], s=35, label="source landmarks after affine")
    ax[0].scatter(points2_yx[:, 1], points2_yx[:, 0], s=35, label="reference landmarks")
    ax[0].set_title("Initial affine source in reference space")
    ax[0].set_aspect("equal")
    ax[0].invert_yaxis()
    ax[0].legend()

    ax[1].imshow(image_for_plot(J), extent=stalign_extent(xJ))
    ax[1].scatter(points2_yx[:, 1], points2_yx[:, 0], s=35, label="reference landmarks")
    ax[1].set_title("Reference H&E")
    ax[1].set_aspect("equal")
    ax[1].invert_yaxis()
    ax[1].legend()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# QC 3: compare final deformed source H&E to reference H&E
def make_stalign_deformed_image_qc_plot(stalign_data, stalign_out, outpath):
    xI = stalign_data["xI"]
    I = stalign_data["I"]
    xJ = stalign_data["xJ"]
    J = stalign_data["J"]

    phiI = STalign.transform_image_source_to_target(
        stalign_out["xv"],
        stalign_out["v"],
        stalign_out["A"],
        xI,
        I,
        xJ,
    )

    fig, ax = plt.subplots(1, 2, figsize=(16, 8))
    ax[0].imshow(image_for_plot(phiI), extent=stalign_extent(xJ))
    ax[0].set_title("Final STaligned source H&E")
    ax[0].set_aspect("equal")
    ax[0].invert_yaxis()

    ax[1].imshow(image_for_plot(J), extent=stalign_extent(xJ))
    ax[1].set_title("Reference H&E")
    ax[1].set_aspect("equal")
    ax[1].invert_yaxis()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# QC 4: plot the nonlinear deformation grid over the final deformed source H&E
def make_stalign_deformation_grid_qc_plot(stalign_data, stalign_out, points1_yx, points2_yx, outpath, grid_levels=20):
    xI = stalign_data["xI"]
    I = stalign_data["I"]
    xJ = stalign_data["xJ"]
    yJ, xJ_axis = xJ

    phii = STalign.build_transform(
        stalign_out["xv"],
        stalign_out["v"],
        stalign_out["A"],
        XJ=xJ,
        direction="b",
    )
    phiI = STalign.transform_image_source_to_target(
        stalign_out["xv"],
        stalign_out["v"],
        stalign_out["A"],
        xI,
        I,
        xJ,
    )
    transformed_landmarks = STalign.transform_points_source_to_target(
        stalign_out["xv"],
        stalign_out["v"],
        stalign_out["A"],
        points1_yx,
    )

    phii = to_numpy(phii)
    transformed_landmarks = to_numpy(transformed_landmarks)

    fig, ax = plt.subplots()
    ax.imshow(image_for_plot(phiI), extent=stalign_extent(xJ), alpha=0.8)
    add_contours(ax, xJ_axis, yJ, phii[..., 0], n_levels=grid_levels)
    add_contours(ax, xJ_axis, yJ, phii[..., 1], n_levels=grid_levels)
    ax.scatter(transformed_landmarks[:, 1], transformed_landmarks[:, 0], s=45, label="source landmarks after STalign")
    ax.scatter(points2_yx[:, 1], points2_yx[:, 0], s=45, label="reference landmarks")
    ax.set_title("STalign deformation grid and landmark fit")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# QC 5: plot aligned Visium source spots over the reference H&E and reference spots
def make_stalign_spots_on_target_qc_plot(src_out, tgt, stalign_data, outpath, title):
    xJ = stalign_data["xJ"]
    J = stalign_data["J"]

    fig, ax = plt.subplots()
    ax.imshow(image_for_plot(J), extent=stalign_extent(xJ), alpha=0.65)
    ax.scatter(tgt["x"], tgt["y"], s=8, alpha=0.35, label="reference spots")
    ax.scatter(src_out["aligned_x"], src_out["aligned_y"], s=8, alpha=0.35, label="source spots after STalign")
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# QC 6: histogram of how far source spots moved after STalign
def make_stalign_displacement_histogram(src_out, outpath):
    if "stalign_total_displacement" not in src_out.columns:
        return
    fig, ax = plt.subplots()
    ax.hist(src_out["stalign_total_displacement"].dropna(), bins=40, alpha=0.55, label="total")
    if "stalign_nonlinear_displacement" in src_out.columns:
        ax.hist(src_out["stalign_nonlinear_displacement"].dropna(), bins=40, alpha=0.55, label="nonlinear after affine")
    ax.set_xlabel("Spot displacement in hires pixels")
    ax.set_ylabel("Number of spots")
    ax.set_title("STalign spot displacement distribution")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# QC 7: WM/matching weights at transformed spot positions --> potentially remove 
# In the notebook these values are used to remove background single cells. For Visium
# spots, they are safer as QC annotation unless you explicitly decide a threshold.
def make_stalign_wm_qc_plots(src_out, spot_outpath, hist_outpath):
    if "stalign_WM_value" not in src_out.columns:
        return

    fig, ax = plt.subplots()
    scatter = ax.scatter(
        src_out["aligned_x"],
        src_out["aligned_y"],
        c=src_out["stalign_WM_value"],
        s=10,
    )
    ax.set_title("STalign WM values at aligned Visium spot positions")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    fig.colorbar(scatter, ax=ax, label="WM value")
    fig.tight_layout()
    fig.savefig(spot_outpath, dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(src_out["stalign_WM_value"].dropna(), bins=30)
    ax.set_xlabel("WM value")
    ax.set_ylabel("Number of spots")
    ax.set_title("Distribution of STalign WM values on Visium spots")
    fig.tight_layout()
    fig.savefig(hist_outpath, dpi=300)
    plt.close(fig)

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
    parser.add_argument("--outdir", default=None)
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
    parser.add_argument("--skip_stalign_qc", action="store_true", help="Skip extra STalign QC plots.")
    parser.add_argument("--stalign_grid_levels", default=20, type=int, help="Number of contour levels for the deformation-grid QC plot.")
    parser.add_argument("--wm_threshold", default=None, type=float, help="Optional WM threshold saved as stalign_WM_pass; the output spots are not filtered.")
    return parser.parse_args() ## validate as in STcompare?

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
        stalign_out, stalign_data = run_stalign_registration(args, points1, points2)

        # applying STalign transform to source spots
        src_out = align_source_spots_stalign(src, stalign_out)

        # add spot-level QC metrics before saving the output table
        src_out = stalign_qc(
            src_out,
            stalign_out,
            stalign_data,
            wm_threshold=args.wm_threshold,
        )

        # landmark fit plot
        make_stalign_landmark_fit_plot(points1, points2, stalign_out, paths["landmark_plot"])

        # extra STalign QC plots adapted from the STalign notebook, but using Visium spots
        # instead of single-cell centroid/rasterized MERFISH data
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
    make_aligned_overlay_plot(
        src_out,
        tgt,
        paths["after_plot"],
        title=f"After alignment: {args.sample_aligned} aligned to {args.sample_reference}",
    )

if __name__ == "__main__":
    main()