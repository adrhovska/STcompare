"""
align_and_add_barcodes.py

Aligns an input Visium sample (aligned) to a reference Visium sample
using STalign LDDMM, then attaches barcodes to the aligned coordinates.

USAGE:
    python align_and_add_barcodes.py \
        --pos1      /path/to/sample_1/spatial/tissue_positions.csv \
        --pos2      /path/to/sample_2/spatial/tissue_positions.csv \
        --scale1    /path/to/sample_1/spatial/scalefactors_json.json \
        --scale2    /path/to/sample_2/spatial/scalefactors_json.json \
        --outdir    /path/to/output \
        --sample_aligned     SampleA \
        --sample_reference     SampleB \
        [--outname  SampleA_aligned_to_SampleB_WITH_BARCODES.csv] \
        [--scale    hires] \
        [--dx       30.0] \
        [--blur     2.0] \
        [--niter    500] \
        [--diffeo_start 100] \
        [--device   cpu]

OUTPUTS:
    <outdir>/<sample_aligned>_aligned_to_<sample_reference>_WITH_BARCODES.csv
        CSV with columns: barcode, x, y
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
import torch
from torch.nn.functional import grid_sample


# =============================================================================
# Core interpolation and transform utilities
# =============================================================================

def _interp(x, I, phii, **kwargs):
    """Interpolate 2D image I at positions phii using torch grid_sample."""
    I = torch.as_tensor(I)
    phii = torch.as_tensor(phii)
    phii = torch.clone(phii)
    for i in range(2):
        phii[i] -= x[i][0]
        phii[i] /= x[i][-1] - x[i][0]
    phii *= 2.0
    phii -= 1.0
    out = grid_sample(
        I[None],
        phii.flip(0).permute((1, 2, 0))[None],
        align_corners=True,
        **kwargs
    )
    return out[0]


def _to_A(L, T):
    """Construct a 3x3 affine matrix from a 2x2 linear part L and translation T."""
    O = torch.tensor([0., 0., 1.], device=L.device, dtype=L.dtype)
    return torch.cat((torch.cat((L, T[:, None]), 1), O[None]))


def _normalize(arr, t_min=0, t_max=1):
    """Linearly normalize an array to [t_min, t_max]."""
    diff_arr = np.max(arr) - np.min(arr)
    if diff_arr == 0:
        return np.zeros_like(arr)
    return ((arr - np.min(arr)) / diff_arr * (t_max - t_min)) + t_min


# =============================================================================
# Rasterization
# =============================================================================

def rasterize_spots(x, y, dx=30.0, blur=1.0, expand=1.1):
    """
    Rasterize Visium spot coordinates into a density image for LDDMM.

    Parameters
    ----------
    x : array-like  — x coordinates of spots
    y : array-like  — y coordinates of spots
    dx : float      — pixel size in same units as x and y
    blur : float    — Gaussian kernel width in pixels
    expand : float  — factor to expand the bounding box

    Returns
    -------
    X_ : numpy array — pixel locations along x axis
    Y_ : numpy array — pixel locations along y axis
    W  : numpy array — rasterized density image, channels on first axis
    """
    blur = [blur] if not isinstance(blur, list) else blur
    blur = np.array(blur)
    maxblur = np.max(blur)

    minx, maxx = np.min(x), np.max(x)
    miny, maxy = np.min(y), np.max(y)
    cx = (minx + maxx) / 2.0
    cy = (miny + maxy) / 2.0
    minx = cx - (maxx - minx) / 2.0 * expand
    maxx = cx + (maxx - minx) / 2.0 * expand
    miny = cy - (maxy - miny) / 2.0 * expand
    maxy = cy + (maxy - miny) / 2.0 * expand

    X_ = np.arange(minx, maxx, dx)
    Y_ = np.arange(miny, maxy, dx)
    X = np.stack(np.meshgrid(X_, Y_))
    W = np.zeros((X.shape[1], X.shape[2], len(blur)))

    r = int(np.ceil(maxblur * 4))

    for x_, y_ in zip(x, y):
        col = int(np.round((x_ - X_[0]) / dx))
        row = int(np.round((y_ - Y_[0]) / dx))
        row0 = max(row - r, 0)
        row1 = min(row + r, W.shape[0] - 1)
        col0 = max(col - r, 0)
        col1 = min(col + r, W.shape[1] - 1)

        k = np.exp(-(
            (X[0][row0:row1+1, col0:col1+1, None] - x_)**2 +
            (X[1][row0:row1+1, col0:col1+1, None] - y_)**2
        ) / (2.0 * (dx * blur * 2)**2))
        k /= np.sum(k, axis=(0, 1), keepdims=True)
        W[row0:row1+1, col0:col1+1, :] += k

    W = np.abs(W).transpose((-1, 0, 1))
    return X_, Y_, W


# =============================================================================
# Position file reading
# =============================================================================

def _has_no_header(path):
    """
    Peek at the first field of a CSV to detect whether a header row is present.
    If the first field looks like a barcode (ACGT characters), there is no header.
    """
    with open(path) as f:
        first = f.readline().split(",")[0].strip()
    return all(c in "ACGTacgt-0123456789" for c in first[:8])


def read_positions(pos_path, scale_path, scale_type="hires", in_tissue_only=True):
    """
    Read a Visium tissue_positions.csv and return a DataFrame with columns:
    barcode, x, y — where x and y are scaled to hires or lowres image space.

    Handles both old 10x format (no header) and new format (with header).

    Parameters
    ----------
    pos_path       : str  — path to tissue_positions.csv
    scale_path     : str  — path to scalefactors_json.json
    scale_type     : str  — 'hires' or 'lowres'
    in_tissue_only : bool — whether to filter to in-tissue spots only

    Returns
    -------
    DataFrame with columns: barcode, x, y
    """
    pos = pd.read_csv(
        pos_path,
        header=None if _has_no_header(pos_path) else 0
    )

    if pos.columns[0] != "barcode":
        pos.columns = [
            "barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres"
        ]

    if "barcode" not in pos.columns:
        pos = pos.reset_index()
        pos.columns = ["barcode"] + list(pos.columns[1:])

    if in_tissue_only and "in_tissue" in pos.columns:
        pos = pos[pos["in_tissue"] == 1].copy()

    with open(scale_path) as f:
        scales = json.load(f)

    sf = scales["tissue_hires_scalef"] if scale_type == "hires" \
        else scales["tissue_lowres_scalef"]

    print(f"  Scale factor ({scale_type}): {sf}")

    pos["x"] = pos["pxl_col_in_fullres"].astype(float) * sf
    pos["y"] = pos["pxl_row_in_fullres"].astype(float) * sf

    return pos[["barcode", "x", "y"]].reset_index(drop=True)


# =============================================================================
# LDDMM alignment
# =============================================================================

def align_samples(pos1, pos2, dx=30.0, blur=2.0,
                  niter=500, diffeo_start=100,
                  a=500.0, epL=2e-8, epT=2e-1, epV=2e3,
                  sigmaM=1.0, sigmaR=5e5,
                  device="cpu", dtype=torch.float64):
    """
    Rasterize both spot sets and run LDDMM to align sample 1 to sample 2.

    Parameters
    ----------
    pos1, pos2    : DataFrames with columns barcode, x, y
    dx            : raster pixel size
    blur          : Gaussian kernel width
    niter         : total number of gradient descent iterations
    diffeo_start  : iteration at which deformable (diffeomorphic) registration begins
    a             : smoothness scale of velocity field
    epL, epT, epV : gradient descent step sizes for linear, translation, velocity
    sigmaM        : image matching weight
    sigmaR        : regularization weight
    device        : torch device string
    dtype         : torch dtype

    Returns
    -------
    A   : torch tensor — 3x3 affine matrix
    v   : torch tensor — velocity field
    xv  : list        — velocity field sample point grids
    xI  : list        — sample 1 pixel grids
    """
    print("  Rasterizing sample 1...")
    xI_, yI_, WI = rasterize_spots(
        pos1["x"].values, pos1["y"].values, dx=dx, blur=blur
    )
    print("  Rasterizing sample 2...")
    xJ_, yJ_, WJ = rasterize_spots(
        pos2["x"].values, pos2["y"].values, dx=dx, blur=blur
    )

    WI = _normalize(WI).astype(np.float64)
    WJ = _normalize(WJ).astype(np.float64)

    xI = [torch.tensor(yI_, dtype=dtype, device=device),
          torch.tensor(xI_, dtype=dtype, device=device)]
    xJ = [torch.tensor(yJ_, dtype=dtype, device=device),
          torch.tensor(xJ_, dtype=dtype, device=device)]

    I = torch.tensor(WI, dtype=dtype, device=device)
    J = torch.tensor(WJ, dtype=dtype, device=device)

    # Affine initialisation
    L = torch.eye(2, dtype=dtype, device=device, requires_grad=True)
    T = torch.zeros(2, dtype=dtype, device=device, requires_grad=True)

    # Velocity field setup
    expand = 2.0
    p = 2.0
    nt = 3
    minv = torch.as_tensor([x[0] for x in xI], dtype=dtype, device=device)
    maxv = torch.as_tensor([x[-1] for x in xI], dtype=dtype, device=device)
    minv = (minv + maxv) * 0.5 - (maxv - minv) * expand * 0.5
    maxv = (minv + maxv) * 0.5 + (maxv - minv) * expand * 0.5
    xv = [torch.arange(m, M, a * 0.5, dtype=dtype, device=device)
          for m, M in zip(minv, maxv)]
    XV = torch.stack(torch.meshgrid(xv, indexing="ij"), -1)
    v = torch.zeros(
        (nt, XV.shape[0], XV.shape[1], XV.shape[2]),
        dtype=dtype, device=device, requires_grad=True
    )

    # Regularization kernel in frequency space
    dv = torch.as_tensor([x[1] - x[0] for x in xv], dtype=dtype, device=device)
    fv = [torch.arange(n, dtype=dtype, device=device) / n / d
          for n, d in zip(XV.shape, dv)]
    FV = torch.stack(torch.meshgrid(fv, indexing="ij"), -1)
    LL = (1.0 + 2.0 * a**2 * torch.sum(
        (1.0 - torch.cos(2.0 * np.pi * FV * dv)) / dv**2, -1
    ))**(p * 2.0)
    K = 1.0 / LL
    DV = torch.prod(dv)

    XJ = torch.stack(torch.meshgrid(*xJ, indexing="ij"), -1)

    print(f"\n  Running LDDMM for {niter} iterations "
          f"(deformable registration starts at iter {diffeo_start})...")

    for it in range(niter):
        A = _to_A(L, T)
        Ai = torch.linalg.inv(A)
        Xs = (Ai[:2, :2] @ XJ[..., None])[..., 0] + Ai[:2, -1]

        for t in range(nt - 1, -1, -1):
            Xs = Xs + _interp(
                xv, -v[t].permute(2, 0, 1), Xs.permute(2, 0, 1)
            ).permute(1, 2, 0) / nt

        AI = _interp(xI, I, Xs.permute(2, 0, 1), padding_mode="border")

        EM = torch.sum((AI - J)**2) / 2.0 / sigmaM**2
        ER = torch.sum(
            torch.sum(
                torch.abs(torch.fft.fftn(v, dim=(1, 2)))**2, dim=(0, -1)
            ) * LL
        ) * DV / 2.0 / v.shape[1] / v.shape[2] / sigmaR**2
        E = EM + ER
        E.backward()

        with torch.no_grad():
            factor = 1.0 / (1.0 + (it >= diffeo_start) * 9)
            L -= epL * factor * L.grad
            T -= epT * factor * T.grad
            L.grad.zero_()
            T.grad.zero_()

            vgrad = torch.fft.ifftn(
                torch.fft.fftn(v.grad, dim=(1, 2)) * K[..., None],
                dim=(1, 2)
            ).real
            if it >= diffeo_start:
                v -= epV * vgrad
            v.grad.zero_()

        if it % 100 == 0 or it == niter - 1:
            print(f"    iter {it:4d}  "
                  f"E={E.item():.4f}  "
                  f"EM={EM.item():.4f}  "
                  f"ER={ER.item():.4f}")

    A = _to_A(L, T).detach()
    print("  Alignment complete.")
    return A, v.detach(), xv, xI


# =============================================================================
# Spot transformation
# =============================================================================

def transform_spots(pos1, A, v, xv, device="cpu", dtype=torch.float64):
    """
    Apply the forward transform (velocity field then affine) to sample 1 spots
    and return a DataFrame with barcodes attached.

    Parameters
    ----------
    pos1   : DataFrame with columns barcode, x, y
    A      : 3x3 affine torch tensor
    v      : velocity field torch tensor
    xv     : list of sample point grids for velocity field
    device : torch device string
    dtype  : torch dtype

    Returns
    -------
    DataFrame with columns: barcode, x, y
    """
    nt = v.shape[0]

    # Points in row-col order (y, x) as expected by _interp
    pts = torch.tensor(
        pos1[["y", "x"]].values,
        dtype=dtype, device=device
    )

    # Integrate velocity field forward
    pts_t = pts.clone()
    for t in range(nt):
        pts_t = pts_t + _interp(
            xv,
            v[t].permute(2, 0, 1),
            pts_t.T[..., None]
        )[..., 0].T / nt

    # Apply affine
    pts_t = (A[:2, :2] @ pts_t.T + A[:2, -1][..., None]).T
    pts_t = pts_t.cpu().numpy()

    return pd.DataFrame({
        "barcode": pos1["barcode"].values,
        "x": pts_t[:, 1],   # col → x
        "y": pts_t[:, 0],   # row → y
    })


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Align a Visium sample (sample 1) to a reference Visium sample (sample 2) "
            "using STalign LDDMM and attach barcodes to the aligned coordinates."
        )
    )
    p.add_argument("--pos1",   required=True,
                   help="tissue_positions.csv for the aligned sample")
    p.add_argument("--pos2",   required=True,
                   help="tissue_positions.csv for the reference sample")
    p.add_argument("--scale1", required=True,
                   help="scalefactors_json.json for the aligned sample")
    p.add_argument("--scale2", required=True,
                   help="scalefactors_json.json for the reference sample")
    p.add_argument("--outdir", required=True,
                   help="Directory to save the output CSV")
    p.add_argument("--name1",  default="sample_1",
                   help="Display name for the aligned sample (default: sample_1)")
    p.add_argument("--name2",  default="sample_2",
                   help="Display name for the reference sample (default: sample_2)")
    p.add_argument("--outname", default=None,
                   help=(
                       "Output filename. "
                       "Defaults to <name1>_aligned_to_<name2>_WITH_BARCODES.csv"
                   ))
    p.add_argument("--scale",  default="hires", choices=["hires", "lowres"],
                   help="Image scale to use for coordinates (default: hires)")
    p.add_argument("--dx",     type=float, default=30.0,
                   help="Raster pixel size in same units as coordinates (default: 30.0)")
    p.add_argument("--blur",   type=float, default=2.0,
                   help="Gaussian blur kernel width in pixels (default: 2.0)")
    p.add_argument("--niter",  type=int,   default=500,
                   help="Number of LDDMM gradient descent iterations (default: 500)")
    p.add_argument("--diffeo_start", type=int, default=100,
                   help="Iteration to begin deformable registration (default: 100)")
    p.add_argument("--device", default="cpu",
                   help="Torch device: cpu or cuda:0 (default: cpu)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    outname = args.outname or (
        f"{args.name1}_aligned_to_{args.name2}_WITH_BARCODES.csv"
    )

    print("\n=== STalign + Barcode Attachment ===")
    print(f"name1        : {args.name1}")
    print(f"name2        : {args.name2}")
    print(f"pos1         : {args.pos1}")
    print(f"pos2         : {args.pos2}")
    print(f"scale        : {args.scale}")
    print(f"outdir       : {args.outdir}")
    print(f"outname      : {outname}")
    print(f"device       : {args.device}")
    print(f"niter        : {args.niter}")
    print(f"diffeo_start : {args.diffeo_start}\n")

    print(f"Reading {args.name1} positions...")
    pos1 = read_positions(args.pos1, args.scale1, scale_type=args.scale)
    print(f"  {args.name1} in-tissue spots: {len(pos1)}")

    print(f"Reading {args.name2} positions...")
    pos2 = read_positions(args.pos2, args.scale2, scale_type=args.scale)
    print(f"  {args.name2} in-tissue spots: {len(pos2)}")

    print(f"\nAligning {args.name1} to {args.name2}...")
    A, v, xv, xI = align_samples(
        pos1, pos2,
        dx=args.dx,
        blur=args.blur,
        niter=args.niter,
        diffeo_start=args.diffeo_start,
        device=args.device,
        dtype=torch.float64
    )

    print(f"\nTransforming {args.name1} spot coordinates...")
    aligned = transform_spots(
        pos1, A, v, xv,
        device=args.device,
        dtype=torch.float64
    )

    print(f"\nAligned spots : {len(aligned)}")
    print("Preview:")
    print(aligned.head())

    outpath = os.path.join(args.outdir, outname)
    aligned.to_csv(outpath, index=False)
    print(f"\nSaved: {outpath}")
    print("=== Done ===\n")


if __name__ == "__main__":
    main()