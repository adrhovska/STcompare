## bc this is Visium spot resolution first the H&E images have to be aligned and this then extrapolated to coordinates 
# load required libraries
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import cv2
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
def infer_hires_image_from_pos(pos_file):
    pos_file = Path(pos_file)
    image_file = pos_file.parent / "tissue_hires_image.png"
    if not image_file.exists():
        raise FileNotFoundError(f"Could not find tissue_hires_image.png next to {pos_file}")
    return image_file

# image read
def read_hires_image(image_file, image_downsample=1):
    image_file = Path(image_file)
    img = plt.imread(image_file)

    # Remove alpha channel if present
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[:, :, :3]

    # If grayscale, add channel dimension
    if img.ndim == 2:
        img = img[:, :, None]

    if image_downsample > 1:
        img = img[::image_downsample, ::image_downsample, :]

    img = img.astype(np.float32)

    # STalign expects channels x y x
    img = np.transpose(img, (2, 0, 1))
    img = STalign.normalize(img)
    Y = np.arange(img.shape[1], dtype=float)
    X = np.arange(img.shape[2], dtype=float)

    return X, Y, img

def read_image_for_feature_matching(image_file, image_downsample=1):
    """
    Read H&E image for OpenCV feature matching.
    Returns grayscale uint8 image in downsampled image coordinates.
    """
    image_file = Path(image_file)
    img = plt.imread(image_file)

    # Remove alpha channel if present
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[:, :, :3]

    # Convert float images from 0-1 to 0-255
    if img.dtype != np.uint8:
        if img.max() <= 1.0:
            img = (img * 255).astype(np.uint8)
        else:
            img = img.astype(np.uint8)

    if image_downsample > 1:
        img = img[::image_downsample, ::image_downsample]

    if img.ndim == 3:
        gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    else:
        gray = img

    # Improve contrast; useful for pale H&E
    clahe = cv2.createCLAHE(
        clipLimit=2.0,
        tileGridSize=(8, 8),
    )
    gray = clahe.apply(gray)

    return gray

#####
def auto_landmarks_orb(
    image1,
    image2,
    image_downsample=1,
    max_landmarks=80,
    min_landmarks=8,
    match_ratio=0.75,
    ransac_thresh=25.0,
    outplot=None,
):
    """
    Automatically detect matching H&E landmarks using ORB + RANSAC.

    Returns:
        pointsI: source points in STalign format [y, x], downsampled image coords
        pointsJ: target/reference points in STalign format [y, x], downsampled image coords
    """
    gray1 = read_image_for_feature_matching(
        image1,
        image_downsample=image_downsample,
    )

    gray2 = read_image_for_feature_matching(
        image2,
        image_downsample=image_downsample,
    )

    orb = cv2.ORB_create(
        nfeatures=8000,
        fastThreshold=5,
    )

    kp1, des1 = orb.detectAndCompute(gray1, None)
    kp2, des2 = orb.detectAndCompute(gray2, None)

    if des1 is None or des2 is None:
        raise RuntimeError("ORB could not find descriptors in one or both images.")

    matcher = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=False)
    raw_matches = matcher.knnMatch(des1, des2, k=2)

    good_matches = []

    for pair in raw_matches:
        if len(pair) < 2:
            continue

        m, n = pair

        if m.distance < match_ratio * n.distance:
            good_matches.append(m)

    if len(good_matches) < min_landmarks:
        raise RuntimeError(
            f"Only {len(good_matches)} raw ORB matches found. "
            f"Need at least {min_landmarks}."
        )

    pts1_xy = np.float32(
        [kp1[m.queryIdx].pt for m in good_matches]
    )

    pts2_xy = np.float32(
        [kp2[m.trainIdx].pt for m in good_matches]
    )

    # RANSAC filters out wrong matches
    affine, inliers = cv2.estimateAffinePartial2D(
        pts1_xy,
        pts2_xy,
        method=cv2.RANSAC,
        ransacReprojThreshold=ransac_thresh,
        maxIters=5000,
        confidence=0.99,
    )

    if affine is None or inliers is None:
        raise RuntimeError("RANSAC failed to estimate an affine transform.")

    inliers = inliers.ravel().astype(bool)

    pts1_xy = pts1_xy[inliers]
    pts2_xy = pts2_xy[inliers]

    if pts1_xy.shape[0] < min_landmarks:
        raise RuntimeError(
            f"Only {pts1_xy.shape[0]} RANSAC inlier matches found. "
            f"Need at least {min_landmarks}."
        )

    # Limit number of landmarks so STalign does not get overloaded
    if pts1_xy.shape[0] > max_landmarks:
        idx = np.linspace(
            0,
            pts1_xy.shape[0] - 1,
            max_landmarks,
        ).astype(int)

        pts1_xy = pts1_xy[idx]
        pts2_xy = pts2_xy[idx]

    # Convert OpenCV x,y to STalign y,x
    pointsI = np.column_stack(
        [
            pts1_xy[:, 1],
            pts1_xy[:, 0],
        ]
    )

    pointsJ = np.column_stack(
        [
            pts2_xy[:, 1],
            pts2_xy[:, 0],
        ]
    )

    print("\nAutomatic landmarks:")
    print(f"raw matches: {len(good_matches)}")
    print(f"RANSAC inliers used: {pointsI.shape[0]}")

    if outplot is not None:
        fig, ax = plt.subplots(figsize=(8, 8))

        ax.imshow(gray1, cmap="gray")
        ax.scatter(
            pointsI[:, 1],
            pointsI[:, 0],
            s=25,
            label="auto source landmarks",
        )
        ax.set_title("Automatic source landmarks")
        ax.invert_yaxis()
        ax.legend()
        fig.tight_layout()
        fig.savefig(outplot, dpi=300)
        plt.close(fig)

    return pointsI, pointsJ, affine


def read_landmark_points(points_file, image_downsample=1):
    """
    Read manual landmark CSV with columns y,x in hires image coordinates.
    Converts to downsampled image coordinates for STalign.
    """
    points_file = Path(points_file)
    df = pd.read_csv(points_file)

    if not {"y", "x"}.issubset(df.columns):
        raise ValueError(f"{points_file} must contain columns named 'y' and 'x'.")

    points = df[["y", "x"]].astype(float).to_numpy()
    points = points / image_downsample

    return points
#####
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
    ax.scatter(tgt["x"], tgt["y"], s=8, alpha=0.45, label="reference spots")
    ax.scatter(src_out["aligned_x_hires"], src_out["aligned_y_hires"], s=8, alpha=0.45, label="source spots aligned")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title("Aligned overlay plot")
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
    parser.add_argument("--coord_scale", default="hires", choices=["um", "fullres", "hires", "lowres"], help=("Coordinate system for alignment. Use 'hires' for H&E image alignment."))
    parser.add_argument("--image1", default=None, help="Source tissue_hires_image.png. If not given, inferred from pos1 folder.")
    parser.add_argument("--image2", default=None, help="Reference tissue_hires_image.png. If not given, inferred from pos2 folder.")
    parser.add_argument("--image_downsample", type=int, default=1, help="Downsample H&E images before alignment. Use 2 or 4 if alignment is slow.")
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
    parser.add_argument("--epL", type=float, default=5e-11)
    parser.add_argument("--epT", type=float, default=5e-4)
    parser.add_argument("--sigmaP", type=float, default=0.2)
    parser.add_argument("--auto_landmarks", action="store_true", help="Automatically detect landmark matches between H&E images using ORB + RANSAC.")
    parser.add_argument("--max_auto_landmarks", type=int, default=80, help="Maximum number of automatic landmarks to pass into STalign.")
    parser.add_argument("--min_auto_landmarks", type=int, default=8, help="Minimum number of automatic landmark matches required.")
    parser.add_argument("--points1", default=None, help="Optional manual source landmarks CSV with columns y,x.")
    parser.add_argument("--points2", default=None, help="Optional manual reference landmarks CSV with columns y,x.")
    parser.add_argument("--auto_match_ratio", type=float, default=0.75, help="Lowe ratio threshold for ORB feature matching.")
    parser.add_argument("--auto_ransac_thresh", type=float, default=25.0, help="RANSAC reprojection threshold in downsampled image pixels.")
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
    # reading H&E images for image-to-image alignment
    if args.image1 is None:
        image1 = infer_hires_image_from_pos(args.pos1)
    else:
        image1 = Path(args.image1)
    if args.image2 is None:
        image2 = infer_hires_image_from_pos(args.pos2)
    else:
        image2 = Path(args.image2)

    XI, YI, I = read_hires_image(image1, image_downsample=args.image_downsample)
    XJ, YJ, J = read_hires_image(image2, image_downsample=args.image_downsample)

    # because the images may be downsampled, spot coordinates must be in the same image scale
    src["x_align"] = src["x"] / args.image_downsample
    src["y_align"] = src["y"] / args.image_downsample
    tgt["x_align"] = tgt["x"] / args.image_downsample
    tgt["y_align"] = tgt["y"] / args.image_downsample
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError(f"Requested device {args.device}, but torch.cuda.is_available() is False.")

    # running LMMD
    print(f"device: {args.device}")
    print(f"niter: {args.niter}")
    print(f"dx: {args.dx}")
    print(f"blur: {args.blur}")

####
    pointsI = None
    pointsJ = None
    L = None
    T = None
    affine_auto = None

    if args.auto_landmarks:
        auto_landmark_plot = outdir / (f"{args.sample_aligned}_auto_landmarks_source.png")
        pointsI, pointsJ, affine_auto = auto_landmarks_orb(
            image1=image1,
            image2=image2,
            image_downsample=args.image_downsample,
            max_landmarks=args.max_auto_landmarks,
            min_landmarks=args.min_auto_landmarks,
            match_ratio=args.auto_match_ratio,
            ransac_thresh=args.auto_ransac_thresh,
            outplot=auto_landmark_plot,
        )

    elif args.points1 is not None and args.points2 is not None:
        pointsI = read_landmark_points(
            args.points1,
            image_downsample=args.image_downsample,
        )

        pointsJ = read_landmark_points(
            args.points2,
            image_downsample=args.image_downsample,
        )

    if pointsI is not None:
        if pointsI.shape != pointsJ.shape:
            raise ValueError("Source and reference landmark arrays must have the same shape.")

        if pointsI.shape[0] < 3:
            raise ValueError("Need at least 3 landmark points.")

        L, T = STalign.L_T_from_points(pointsI, pointsJ)

        print("\nUsing landmark initialization")
        print(f"Number of landmarks: {pointsI.shape[0]}")
####

    params = {"niter": args.niter, "device": args.device, "diffeo_start": args.diffeo_start, "sigmaM": args.sigmaM, "sigmaB": args.sigmaB, "sigmaA": args.sigmaA,
            "epV": args.epV, "epL": args.epL, "epT": args.epT}

    if pointsI is not None:
        params.update(
            {"pointsI": pointsI, "pointsJ": pointsJ, "L": L, "T": T, "sigmaP": args.sigmaP}
        )

    out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **params)
    ####
    print("\nAffine matrix A:")
    print(tensor_to_numpy(out["A"]))
    ####
    A = out["A"]
    v = out["v"]
    xv = out["xv"]
    # dource spot transformation
    # transform source spot coordinates using the image-derived transform
    # STalign wants points as [y, x], not [x, y]
    source_points_yx = np.stack([src["y_align"].to_numpy(), src["x_align"].to_numpy()], axis=1)
    transformed_yx = STalign.transform_points_source_to_target(xv, v, A, source_points_yx)
    transformed_yx = tensor_to_numpy(transformed_yx)
    src_out = src.copy()
    # original coordinates in hires image space
    src_out["original_x_hires"] = src_out["x"]
    src_out["original_y_hires"] = src_out["y"]
    # aligned coordinates returned in downsampled image space
    src_out["aligned_y_downsampled"] = transformed_yx[:, 0]
    src_out["aligned_x_downsampled"] = transformed_yx[:, 1]
    # convert aligned coordinates back to reference hires image space
    src_out["aligned_y_hires"] = src_out["aligned_y_downsampled"] * args.image_downsample
    src_out["aligned_x_hires"] = src_out["aligned_x_downsampled"] * args.image_downsample
    # x and y become aligned coordinates in the reference hires image system
    src_out["x"] = src_out["aligned_x_hires"]
    src_out["y"] = src_out["aligned_y_hires"]
    # overwrite x and y as aligned coords
    src_out["x"] = src_out["aligned_x_hires"]
    src_out["y"] = src_out["aligned_y_hires"]

    preferred_cols = [
        "barcode",
        "x",
        "y",
        "original_x_hires",
        "original_y_hires",
        "aligned_x_hires",
        "aligned_y_hires",
        "aligned_x_downsampled",
        "aligned_y_downsampled",
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
