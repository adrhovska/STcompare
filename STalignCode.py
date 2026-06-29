"""
align_and_add_barcodes.py

Aligns an input Visium sample (aligned) to a reference Visium sample
using STalign LDDMM, then attaches barcodes to the aligned coordinates. --> or check barcode presence, fix ?
Source: Arrays of x and y positions of cells from single-cell resolution ST data
Target: Registered H&E image from spot-resolution ST data

USAGE:
    python STalign.py \
        --pos1      /path/to/sample_1/spatial/tissue_positions.csv \
        --pos2      /path/to/sample_2/spatial/tissue_positions.csv \
        --scale1    /path/to/sample_1/spatial/scalefactors_json.json \
        --scale2    /path/to/sample_2/spatial/scalefactors_json.json \
        --outdir    /path/to/output \
        --sample_aligned     SampleA \
        --sample_reference     SampleB \
        [--outname  SampleA_aligned_to_SampleB_barcodes.csv] \
        [--scale    hires] \
        [--dx       30.0] \
        [--blur     2.0] \
        [--niter    500] \
        [--diffeo_start 100] \
        [--device   cpu]

OUTPUTS:
    <outdir>/<sample_aligned>_aligned_to_<sample_reference>_barcodes.csv
        CSV with columns: barcode, x, y
"""
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
plt.rcParams["figure.figsize"] = (12,10)
# make export directory
make_dir= Path("/Users/adrhovska/Desktop/STdata/STcompare_code")
# expected positions
expected_cols = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]
# reading Visium spot data
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

# reading scale factors 
def read_scalefactors(scale_file):
    scale_file = Path(scale_file)
    with open(scale_file, "r") as f:
        scalefactors = json.load(f)
    return scalefactors
# pick and return multiplier for full-res coordinates 
def coordinate_factor(scalefactors, coord_scale, spot_diameter_um):
    if coord_scale == "fullres":
        return 1.0
    if coord_scale == "hires":
        return float(scalefactors["tissue_hires_scalef"])
    if coord_scale == "lowres":
        return float(scalefactors["tissue_lowres_scalef"])
    if coord_scale == "um":
        if "spot_diameter_fullres" not in scalefactors:
            raise KeyError("Could not find 'spot_diameter_fullres' in scalefactors_json.json. ")
        return float(spot_diameter_um) / float(scalefactors["spot_diameter_fullres"])
    raise ValueError(f"Unknown coord_scale: {coord_scale}")
# adding coordinates in the chosen system 
def add_xy_coordinates(df, scalefactors, coord_scale="um", spot_diameter_um=55.0):
    df = df.copy()
    factor = coordinate_factor(
        scalefactors=scalefactors,
        coord_scale=coord_scale,
        spot_diameter_um=spot_diameter_um,
    )
    df["x"] = df["pxl_col_in_fullres"].astype(float) * factor
    df["y"] = df["pxl_row_in_fullres"].astype(float) * factor
    return df
# filter tissue spots in the coord system 
def filter_spots(df):
    df = df.copy()
    df = df[df["in_tissue"].astype(int) == 1].copy()
    df = df.reset_index(drop=True)
    return df
# validation 
def validate_spots(df, sample_name):
    if df.empty:
        raise ValueError(f"{sample_name}: no in-tissue spots found")
    if df["barcode"].duplicated().any():
        duplicated = df.loc[df["barcode"].duplicated(), "barcode"].head().tolist()
        raise ValueError(
            f"{sample_name}: duplicate barcodes found"
        )
    if not np.isfinite(df[["x", "y"]].to_numpy()).all():
        raise ValueError(f"{sample_name}: non-finite x/y coordinates found")
    print(f"{sample_name}: {df.shape[0]} in-tissue spots")
    print(
        f"{sample_name}: x range {df['x'].min():.2f} to {df['x'].max():.2f}, "
        f"y range {df['y'].min():.2f} to {df['y'].max():.2f}"
    )

def tensor_to_numpy(x):
    if torch.is_tensor(x):
        if x.is_cuda:
            x = x.cpu()
        x = x.detach().numpy()
    return np.asarray(x)


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
    ax.scatter(tgt["x"], tgt["y"], s=8, alpha=0.45, label="reference")
    ax.scatter(src_out["aligned_x"], src_out["aligned_y"], s=8, alpha=0.45, label="source aligned")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

def make_raster_plot(XI, YI, I, XJ, YJ, J, outpath):
    extentI = STalign.extent_from_x((YI, XI))
    extentJ = STalign.extent_from_x((YJ, XJ))
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].imshow(I[0], extent=extentI)
    ax[0].invert_yaxis()
    ax[0].set_title("Source raster")
    ax[1].imshow(J[0], extent=extentJ)
    ax[1].invert_yaxis()
    ax[1].set_title("Reference raster")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# saving STalign transform outputs in numpy format
def save_transform_npz(out, outpath):
    to_save = {}
    for key in ["A", "v", "xv"]:
        if key in out:
            to_save[key] = tensor_to_numpy(out[key])
    np.savez_compressed(outpath, **to_save)

# parsing arguments 
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pos1", required=True, help="Source tissue_positions.csv")
    parser.add_argument("--pos2", required=True, help="Reference tissue_positions.csv")
    parser.add_argument("--scale1", required=True, help="Source scalefactors_json.json")
    parser.add_argument("--scale2", required=True, help="Reference scalefactors_json.json")
    parser.add_argument("--project_dir", default=str(make_dir), help="Main project directory where STalign output folders will be created.")
    parser.add_argument("--outdir", default=None, help="Optional custom output directory.")
    parser.add_argument("--sample_aligned", required=True)
    parser.add_argument("--sample_reference", required=True)
    parser.add_argument("--outname", default=None)
    parser.add_argument("--coord_scale", default="um", choices=["um", "fullres", "hires", "lowres"], help=(
            "Coordinate system for alignment. "
            "Use 'um' for cross-sample Visium spot alignment unless you have a reason not to."
        ),
    )
    parser.add_argument("--spot_diameter_um", type=float, default=55.0, help="Standard Visium spot diameter is usually 55 um.",)
    parser.add_argument("--dx", type=float, default=50.0)
    parser.add_argument("--blur", type=float, default=1.5)
    parser.add_argument("--niter", type=int, default=1000)
    parser.add_argument("--diffeo_start", type=int, default=100)
    parser.add_argument("--device", default="cpu", help="Use 'cpu' or 'cuda:0'. STalign can be temperamental")
    parser.add_argument("--sigmaM", type=float, default=0.2)
    parser.add_argument("--sigmaB", type=float, default=2.0)
    parser.add_argument("--sigmaA", type=float, default=5.0)
    parser.add_argument("--epV", type=float, default=50.0)

    return parser.parse_args()

def main():
    args = parse_args()
    project_dir = Path(args.project_dir)
    if args.outdir is None:
        outdir = (project_dir
            / "STalign_outputs"
            / f"{args.sample_aligned}_aligned_to_{args.sample_reference}"
        )
    else: outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.outname is None:
        outname = (
            f"{args.sample_aligned}_aligned_to_"
            f"{args.sample_reference}_barcodes.csv"
        )
    else:
        outname = args.outname
    output_csv = outdir / outname

    # reading source sample
    src_raw = read_spots(args.pos1)
    src_sf = read_scalefactors(args.scale1)

    # reading reference sample
    tgt_raw = read_spots(args.pos2)
    tgt_sf = read_scalefactors(args.scale2)

    # adding coords 
    src = add_xy_coordinates(src_raw, src_sf, coord_scale=args.coord_scale, spot_diameter_um=args.spot_diameter_um)
    tgt = add_xy_coordinates(tgt_raw, tgt_sf, coord_scale=args.coord_scale, spot_diameter_um=args.spot_diameter_um)
    src = filter_spots(src)
    tgt = filter_spots(tgt)
    validate_spots(src, args.sample_aligned)
    validate_spots(tgt, args.sample_reference)
    before_plot = outdir / (f"{args.sample_aligned}_vs_{args.sample_reference}_before_alignment.png")
    # overlay plot generation
    make_overlay_plot(src, tgt, before_plot, title=f"Before alignment: {args.sample_aligned} vs {args.sample_reference}",)
    # rasterising spot coordinates
    xI = src["x"].to_numpy()
    yI = src["y"].to_numpy()
    xJ = tgt["x"].to_numpy()
    yJ = tgt["y"].to_numpy()
    XI, YI, I, figI = STalign.rasterize(xI, yI, dx=args.dx, blur=args.blur)
    XJ, YJ, J, figJ = STalign.rasterize(xJ, yJ, dx=args.dx, blur=args.blur)
    raster_plot = outdir / (f"{args.sample_aligned}_to_{args.sample_reference}_rasters.png")
    make_raster_plot(XI, YI, I, XJ, YJ, J, raster_plot)
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError(f"Requested device {args.device}, but torch.cuda.is_available() is False.")

    # running LMMD
    print(f"device: {args.device}")
    print(f"niter: {args.niter}")
    print(f"dx: {args.dx}")
    print(f"blur: {args.blur}")

    params = {"niter": args.niter, "device": args.device, "diffeo_start": args.diffeo_start, "sigmaM": args.sigmaM, "sigmaB": args.sigmaB,
            "sigmaA": args.sigmaA, "epV": args.epV}
    out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **params)
    A = out["A"]
    v = out["v"]
    xv = out["xv"]
    # dource spot transformation
    source_points_yx = np.stack([src["y"].to_numpy(), src["x"].to_numpy()], axis=1)
    transformed_yx = STalign.transform_points_source_to_target(xv, v, A, source_points_yx)
    transformed_yx = tensor_to_numpy(transformed_yx)
    src_out = src.copy()
    src_out["original_x"] = src_out["x"]
    src_out["original_y"] = src_out["y"]
    src_out["aligned_y"] = transformed_yx[:, 0]
    src_out["aligned_x"] = transformed_yx[:, 1]
    # overwrite x and y as aligned coords
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
    # making afteralignmen plot
    after_plot = outdir / (f"{args.sample_aligned}_aligned_to_{args.sample_reference}_after_alignment.png")
    make_aligned_overlay_plot(src_out, tgt, after_plot, title=f"After alignment: {args.sample_aligned} aligned to {args.sample_reference}",)
    transform_file = outdir / (f"{args.sample_aligned}_to_{args.sample_reference}_STalign_transform.npz")

    save_transform_npz(out, transform_file)

if __name__ == "__main__":
    main()