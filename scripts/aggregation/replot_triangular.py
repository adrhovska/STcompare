"""
Re-plots the aggregated similarity/correlation matrices as triangular
heatmaps (half the square hidden), with the solid right-angle corner
positioned at the bottom-right -- matching the lab's reference heatmap
layout (all_organoids_spatial_similarity_heatmap.pdf), where:
  - rows (y-axis, top->bottom) are in the standard descending sample order
  - columns (x-axis, left->right) are the SAME order reversed
  - a cell is shown only if (row_position + col_position) >= n-1
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

IN_DIR = sys.argv[1]
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else IN_DIR
os.makedirs(OUT_DIR, exist_ok=True)


def triangular_display(mat):
    samples = list(mat.index)  # already in descending reference order
    n = len(samples)
    row_order = samples
    col_order = samples[::-1]
    disp = mat.loc[row_order, col_order]
    mask = np.fromfunction(lambda i, j: (i + j) < (n - 1), (n, n))
    return disp, mask


def plot_triangle(mat, title, out_pdf, cmap, vmin, vmax, center, label):
    disp, mask = triangular_display(mat)
    plt.figure(figsize=(11, 10))
    ax = sns.heatmap(
        disp, mask=mask, vmin=vmin, vmax=vmax, center=center, cmap=cmap,
        square=True, annot=True, fmt=".2f", annot_kws={"size": 6},
        cbar_kws={"label": label}, linewidths=0.4, linecolor="white",
    )
    ax.set_facecolor("white")
    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()
    print("wrote", out_pdf)


# 1. overall similarity (sequential, 0 to max off-diagonal)
overall = pd.read_csv(os.path.join(IN_DIR, "All_Pairwise_Similarity_Matrix.csv"), index_col=0)
off = overall.values[~np.eye(len(overall), dtype=bool)]
vmax = np.nanmax(off)
plot_triangle(
    overall, "Overall spatial similarity (all genes)",
    os.path.join(OUT_DIR, "spatial_similarity_heatmap_triangle.pdf"),
    cmap="magma", vmin=0, vmax=vmax, center=None, label="Similarity",
)

# 2. cluster-level correlation matrices (diverging, -1 to 1, centered at 0)
for cname in ["progenitor_genes", "maturation_genes", "patterning_genes"]:
    f = os.path.join(IN_DIR, f"Cluster_{cname}_Correlation_Matrix.csv")
    mat = pd.read_csv(f, index_col=0)
    plot_triangle(
        mat, f"Cluster-level spatial correlation: {cname}",
        os.path.join(OUT_DIR, f"Cluster_{cname}_correlation_heatmap_triangle.pdf"),
        cmap="RdBu_r", vmin=-1, vmax=1, center=0, label="Correlation coefficient",
    )
