# QC plotting functions for STalignCode.py

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import torch
from STalign import STalign

# make plots bigger uniformly
plt.rcParams["figure.figsize"] = (12, 10)

#// function converting torch tensors (used internally by STalign) to numpy (used in table code and plots)
# due to STalign using torch tensors but the rest of the code using numpy arrays
def to_numpy(x):
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)

#// function preparing images for matplotlib plotting 
# converts to numpy, makes channel-last, squeezes single-channel images and rescales on 0 to 1 range
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

#// function getting the extent of the image for plotting in matplotlib
def stalign_extent(x_axes):
    try:
        return STalign.extent_from_x(tuple(x_axes))
    except Exception:
        y_axis, x_axis = x_axes
        return [x_axis[0], x_axis[-1], y_axis[-1], y_axis[0]]

#// function converting STalign image and axes to matplotlib format for plotting
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


#// function applying affine transformation to x and y coordinates (landmarks and source spots)
def apply_affine(x, y, affine):
    xy_hom = np.column_stack([x, y, np.ones(len(x))])
    transformed_xy = xy_hom @ affine.T 
    return transformed_xy[:, 0], transformed_xy[:, 1]


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

## STalign QC plots 

#// QC1: function making before-alignment plot
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

#// QC2: function making after-alignment plot
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

#// QC3: function to generate landmark fit plot (to show how the reference landmarks align to the source landmarks)
def make_stalign_landmark_fit_plot(points1_yx, points2_yx, stalign_out, outpath):
    transformed_yx = STalign.transform_points_source_to_target(
        stalign_out["xv"],
        stalign_out["v"],
        stalign_out["A"],
        points1_yx,
    )

# making numpy arrays for QC (now in x an y order for plotting)
    transformed_yx = to_numpy(transformed_yx)
    reference_x = points2_yx[:, 1]
    reference_y = points2_yx[:, 0]
    aligned_x = transformed_yx[:, 1]
    aligned_y = transformed_yx[:, 0]

# calculate residuals, the distance between transformed and reference point 
    residuals = np.linalg.norm(transformed_yx - points2_yx, axis=1)

# make the figure (including source and referece landmarks with pair lines/residuals)
    fig, ax = plt.subplots()
    ax.scatter(reference_x, reference_y, s=45, label="reference landmarks")
    ax.scatter(aligned_x, aligned_y, s=45, label="source landmarks after alignment")
    for x1, y1, x2, y2 in zip(aligned_x, aligned_y, reference_x, reference_y):
        ax.plot([x1, x2], [y1, y2], linewidth=0.8, alpha=0.7)
    ax.set_aspect("equal")
    ax.invert_yaxis() # inverted because top left corner is the start of axis otherwise (y increases down, has to be inverted)
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

#// QC4: function for generating LDDMM diagnostic plot 
# pulls out prepared images and their coordinate systems as well as alignment outputs
def make_stalign_lddmm_diagnostic_plot(stalign_data, stalign_out, outpath):
    xI = stalign_data["xI"]
    I = stalign_data["I"]
    xJ = stalign_data["xJ"]
    J = stalign_data["J"]
    xv = stalign_out["xv"]
    v = stalign_out["v"]
    A = stalign_out["A"]
    WM = stalign_out.get("WM", None)

    # phiI is the generated final transformed source image (affine, nonlinear LDDMM both applied as well as placing into reference coordinate system)
    phiI = STalign.transform_image_source_to_target(xv, v, A, xI, I, xJ)

    # phii is the coordinate deformation map on the target grid (backward warping used)
    phii = STalign.build_transform(xv, v, A, XJ=xJ, direction="b")

    # using previously def function to prepare for matplotlib plotting
    J_plot = image_for_plot(J)
    phiI_plot = image_for_plot(phiI)

    # make panels use same reference coordinate system 
    # xJ_axis is used because xJ is used for the whole coordinate system 
    yJ, xJ_axis = xJ
    extentJ = stalign_extent(xJ)

    # residual error calculation and plotting (J is target image, phiI is transformed source image)
    # can handle both greyscale and colored images 
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

    # calculation and plotting of approximate deformation magnitude (how much each pixel has moved)
        yy, xx = np.meshgrid(yJ, xJ_axis, indexing="ij")
        deformation_mag = np.sqrt(
            (phii_np[..., 0] - yy) ** 2 +
            (phii_np[..., 1] - xx) ** 2
        )
        im = ax[1, 0].imshow(deformation_mag, extent=extentJ)
        ax[1, 0].set_title("Deformation magnitude")
        fig.colorbar(im, ax=ax[1, 0], fraction=0.046, pad=0.04)

    # plotting of WM diagnostic map
    if WM is not None:
        WM_np = to_numpy(WM)
        im = ax[1, 1].imshow(WM_np, extent=extentJ)
        ax[1, 1].set_title("STalign weights / WM")
        fig.colorbar(im, ax=ax[1, 1], fraction=0.046, pad=0.04)
    else:
        ax[1, 1].axis("off")
        ax[1, 1].set_title("No WM found")

    # plotting transformed sourcea and deformation grid 
    ax[1, 2].imshow(phiI_plot, extent=extentJ)
    if phii_np.ndim >= 3:
        add_contours(ax[1, 2], xJ_axis, yJ, phii_np[..., 0], n_levels=20)
        add_contours(ax[1, 2], xJ_axis, yJ, phii_np[..., 1], n_levels=20)
    ax[1, 2].set_title("Transformed source with deformation grid")

    # formatting all QC plot panels 
    for a in ax.ravel():
        a.set_aspect("equal")
        a.invert_yaxis()
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

# QC5: visualize the exact H&E images (ovrelapping) and landmarks given to STalign
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

# QC6: affine alignment before LDDMM
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

# QC7: compare final deformed source H&E to reference H&E
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

# QC8: plot the nonlinear deformation grid over the final deformed source H&E
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

# QC9: plot aligned Visium source spots over the reference H&E and reference spots
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

# QC10: histogram of how far source spots moved after STalign
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

# QC11: WM/matching weights at transformed spot positions --> potentially remove 
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

#// QCsum: function making one summary panel with all generated QC plots 
def make_qc_summary_panel(paths, alignment_method, outpath, title=None):
    qc_order = [
        ("QC1", "Before alignment", "before_plot", ["affine", "stalign"]),
        ("QC2", "Input H&E + landmarks", "stalign_input_qc_plot", ["stalign"]),
        ("QC3", "Initial affine", "stalign_initial_affine_qc_plot", ["stalign"]),
        ("QC4", "Landmark fit", "landmark_plot", ["affine", "stalign"]),
        ("QC5", "Deformed source vs reference", "stalign_deformed_image_qc_plot", ["stalign"]),
        ("QC6", "Deformation grid", "stalign_deformation_grid_qc_plot", ["stalign"]),
        ("QC7", "Aligned spots on target", "stalign_spots_on_target_qc_plot", ["stalign"]),
        ("QC8", "Spot displacement histogram", "stalign_displacement_hist_plot", ["stalign"]),
        ("QC9", "LDDMM diagnostic", "stalign_lddmm_diagnostic_plot", ["stalign"]),
        ("QC10", "After alignment", "after_plot", ["affine", "stalign"]),
        ("QC11", "WM values on spots", "stalign_wm_spot_plot", ["stalign"]),
        ("QC12", "WM value histogram", "stalign_wm_hist_plot", ["stalign"]),
    ]

    existing_qc = []
    for qc_number, qc_name, path_key, methods in qc_order:
        if alignment_method not in methods:
            continue

        path = paths.get(path_key, None)
        if path is not None and Path(path).is_file():
            existing_qc.append((qc_number, qc_name, Path(path)))

    if not existing_qc:
        print("No QC plots found for summary panel")
        return

    n_plots = len(existing_qc)
    ncols = 3
    nrows = int(np.ceil(n_plots / ncols))

    fig, ax = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))

    ax = np.asarray(ax).ravel()

    for i, (qc_number, qc_name, path) in enumerate(existing_qc):
        img = plt.imread(path)
        ax[i].imshow(img)
        ax[i].axis("off")
        ax[i].set_title(f"{qc_number}: {qc_name}", fontsize=11)

    for j in range(n_plots, len(ax)):
        ax[j].axis("off")

    if title is not None:
        fig.suptitle(title, fontsize=16)

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)