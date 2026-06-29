# Manual-landmark affine alignment for Visium spot coordinates
# The H&E images are aligned using manually selected landmark pairs,
# then the fitted affine transform is applied to the Visium spot coordinates.

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (12, 10)

make_dir = Path("/Users/adrhovska/Desktop/STdata/STcompare_code")

expected_cols = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]


def read_spots(pos_file):
    pos_file = Path(pos_file)
    df = pd.read_csv(pos_file)

    if not set(expected_cols).issubset(df.columns):
        df = pd.read_csv(
            pos_file,
            header=None,
            names=expected_cols,
        )

    df["barcode"] = df["barcode"].astype(str)

    numeric_cols = [
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def read_scalefactors(scale_file):
    scale_file = Path(scale_file)

    with open(scale_file, "r") as f:
        scalefactors = json.load(f)

    return scalefactors


def coordinate_factor(scalefactors, coord_scale, spot_diameter_um):
    if coord_scale == "fullres":
        return 1.0

    if coord_scale == "hires":
        return float(scalefactors["tissue_hires_scalef"])

    if coord_scale == "lowres":
        return float(scalefactors["tissue_lowres_scalef"])

    if coord_scale == "um":
        if "spot_diameter_fullres" not in scalefactors:
            raise KeyError(
                "Could not find 'spot_diameter_fullres' in scalefactors_json.json."
            )

        return float(spot_diameter_um) / float(scalefactors["spot_diameter_fullres"])

    raise ValueError(f"Unknown coord_scale: {coord_scale}")


def add_xy_coordinates(df, scalefactors, coord_scale="hires", spot_diameter_um=55.0):
    df = df.copy()

    factor = coordinate_factor(
        scalefactors=scalefactors,
        coord_scale=coord_scale,
        spot_diameter_um=spot_diameter_um,
    )

    df["x"] = df["pxl_col_in_fullres"].astype(float) * factor
    df["y"] = df["pxl_row_in_fullres"].astype(float) * factor

    return df


def filter_spots(df):
    df = df.copy()
    df = df[df["in_tissue"].astype(int) == 1].copy()
    df = df.reset_index(drop=True)

    return df


def validate_spots(df, sample_name):
    if df.empty:
        raise ValueError(f"{sample_name}: no in-tissue spots found")

    if df["barcode"].duplicated().any():
        raise ValueError(f"{sample_name}: duplicate barcodes found")

    if not np.isfinite(df[["x", "y"]].to_numpy()).all():
        raise ValueError(f"{sample_name}: non-finite x/y coordinates found")

    print(f"{sample_name}: {df.shape[0]} in-tissue spots")
    print(
        f"{sample_name}: x range {df['x'].min():.2f} to {df['x'].max():.2f}, "
        f"y range {df['y'].min():.2f} to {df['y'].max():.2f}"
    )


def read_landmark_points(points_file):
    """
    Read manual landmark CSV with columns y,x.
    These coordinates must be in the same coordinate system as --coord_scale.

    For normal H&E landmark picking from tissue_hires_image.png, use:
        --coord_scale hires
    """
    points_file = Path(points_file)
    df = pd.read_csv(points_file)

    if not {"y", "x"}.issubset(df.columns):
        raise ValueError(f"{points_file} must contain columns named 'y' and 'x'.")

    points = df[["y", "x"]].astype(float).to_numpy()

    if not np.isfinite(points).all():
        raise ValueError(f"{points_file} contains non-finite landmark coordinates.")

    return points


def fit_affine_from_landmarks(points1_yx, points2_yx):
    """
    Fit affine transform from source landmarks to reference landmarks.

    Input landmark arrays are [y, x].
    Returned affine maps [x, y, 1] in source coordinates to [x, y] in reference coordinates.
    """
    if points1_yx.shape != points2_yx.shape:
        raise ValueError("Source and reference landmark arrays must have the same shape.")

    if points1_yx.shape[0] < 3:
        raise ValueError("Need at least 3 landmark pairs to fit an affine transform.")

    source_xy = np.column_stack([points1_yx[:, 1], points1_yx[:, 0]])
    reference_xy = np.column_stack([points2_yx[:, 1], points2_yx[:, 0]])

    design = np.column_stack([source_xy, np.ones(source_xy.shape[0])])

    coef_x, *_ = np.linalg.lstsq(design, reference_xy[:, 0], rcond=None)
    coef_y, *_ = np.linalg.lstsq(design, reference_xy[:, 1], rcond=None)

    affine = np.vstack([coef_x, coef_y])

    predicted_xy = design @ affine.T
    residuals = np.linalg.norm(predicted_xy - reference_xy, axis=1)

    return affine, residuals


def apply_affine_to_xy(x, y, affine):
    xy_hom = np.column_stack([x, y, np.ones(len(x))])
    transformed_xy = xy_hom @ affine.T

    return transformed_xy[:, 0], transformed_xy[:, 1]


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


def make_landmark_fit_plot(points1_yx, points2_yx, affine, outpath):
    source_x = points1_yx[:, 1]
    source_y = points1_yx[:, 0]
    reference_x = points2_yx[:, 1]
    reference_y = points2_yx[:, 0]

    aligned_x, aligned_y = apply_affine_to_xy(source_x, source_y, affine)

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


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--pos1", required=True, help="Source tissue_positions.csv")
    parser.add_argument("--pos2", required=True, help="Reference tissue_positions.csv")
    parser.add_argument("--scale1", required=True, help="Source scalefactors_json.json")
    parser.add_argument("--scale2", required=True, help="Reference scalefactors_json.json")
    parser.add_argument("--points1", required=True, help="Source manual landmarks CSV with columns y,x")
    parser.add_argument("--points2", required=True, help="Reference manual landmarks CSV with columns y,x")

    parser.add_argument(
        "--project_dir",
        default=str(make_dir),
        help="Main project directory where output folders will be created.",
    )
    parser.add_argument("--outdir", default=None, help="Optional custom output directory.")
    parser.add_argument("--sample_aligned", required=True)
    parser.add_argument("--sample_reference", required=True)
    parser.add_argument("--outname", default=None)
    parser.add_argument(
        "--coord_scale",
        default="hires",
        choices=["um", "fullres", "hires", "lowres"],
        help="Coordinate system for spot coordinates and landmark coordinates.",
    )
    parser.add_argument(
        "--spot_diameter_um",
        type=float,
        default=55.0,
        help="Standard Visium spot diameter is usually 55 um.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    project_dir = Path(args.project_dir)

    if args.outdir is None:
        outdir = (
            project_dir
            / "STalign_outputs"
            / f"{args.sample_aligned}_aligned_to_{args.sample_reference}"
        )
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

    output_csv = outdir / outname

    src_raw = read_spots(args.pos1)
    src_sf = read_scalefactors(args.scale1)

    tgt_raw = read_spots(args.pos2)
    tgt_sf = read_scalefactors(args.scale2)

    src = add_xy_coordinates(
        src_raw,
        src_sf,
        coord_scale=args.coord_scale,
        spot_diameter_um=args.spot_diameter_um,
    )
    tgt = add_xy_coordinates(
        tgt_raw,
        tgt_sf,
        coord_scale=args.coord_scale,
        spot_diameter_um=args.spot_diameter_um,
    )

    src = filter_spots(src)
    tgt = filter_spots(tgt)

    validate_spots(src, args.sample_aligned)
    validate_spots(tgt, args.sample_reference)

    points1 = read_landmark_points(args.points1)
    points2 = read_landmark_points(args.points2)

    affine, residuals = fit_affine_from_landmarks(points1, points2)

    print("\nManual landmark affine matrix:")
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

    before_plot = outdir / (
        f"{args.sample_aligned}_vs_{args.sample_reference}_before_alignment.png"
    )
    make_overlay_plot(
        src,
        tgt,
        before_plot,
        title=f"Before alignment: {args.sample_aligned} vs {args.sample_reference}",
    )

    aligned_x, aligned_y = apply_affine_to_xy(src["x"].to_numpy(), src["y"].to_numpy(), affine)

    src_out = src.copy()
    src_out["original_x"] = src_out["x"]
    src_out["original_y"] = src_out["y"]
    src_out["aligned_x"] = aligned_x
    src_out["aligned_y"] = aligned_y
    src_out["x"] = src_out["aligned_x"]
    src_out["y"] = src_out["aligned_y"]

    preferred_cols = [
        "barcode",
        "x",
        "y",
        "original_x",
        "original_y",
        "aligned_x",
        "aligned_y",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]

    remaining_cols = [c for c in src_out.columns if c not in preferred_cols]
    src_out = src_out[preferred_cols + remaining_cols]
    src_out.to_csv(output_csv, index=False)

    after_plot = outdir / (
        f"{args.sample_aligned}_aligned_to_{args.sample_reference}_after_alignment.png"
    )
    make_aligned_overlay_plot(
        src_out,
        tgt,
        after_plot,
        title=f"After alignment: {args.sample_aligned} aligned to {args.sample_reference}",
    )

    landmark_plot = outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_manual_landmark_fit.png"
    )
    make_landmark_fit_plot(points1, points2, affine, landmark_plot)

    transform_file = outdir / (
        f"{args.sample_aligned}_to_{args.sample_reference}_manual_affine_transform.npz"
    )
    np.savez_compressed(
        transform_file,
        affine=affine,
        landmark_residuals=residuals,
        points1=points1,
        points2=points2,
        coord_scale=args.coord_scale,
    )

    print("\nDone with manual-landmark affine alignment.")
    print(f"Aligned coordinates: {output_csv}")
    print(f"Before plot: {before_plot}")
    print(f"After plot: {after_plot}")
    print(f"Landmark fit plot: {landmark_plot}")
    print(f"Transform file: {transform_file}")


if __name__ == "__main__":
    main()
