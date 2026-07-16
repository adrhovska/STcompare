"""
Reconstructs the lab's proprietary similarity matrix from values transcribed
(cell-by-cell, cross-checked against multiple high-resolution crops of the
vector PDF) out of all_organoids_spatial_similarity_heatmap.pdf.

One cell -- (Donor2_day40_dif1_org2, Donor2_day120_dif2_org2) -- renders as
pure black with no legible text in the source PDF itself (not a
transcription/resolution issue: pixel-level inspection shows the tile fill
is exactly RGB(0,0,0) with no anti-aliased glyph edges anywhere in the
cell), so it is recorded as NaN and excluded from downstream analysis.
"""
import numpy as np
import pandas as pd

samples = [
    "Donor2_day120_dif2_org2", "Donor2_day120_dif2_org1", "Donor2_day120_dif1_org2",
    "Donor2_day120_dif1_org1", "Donor1_day120_dif2_org2", "Donor1_day120_dif2_org1",
    "Donor1_day120_dif1_org2", "Donor1_day120_dif1_org1", "Donor2_day70_dif2_org2",
    "Donor2_day70_dif2_org1", "Donor2_day70_dif1_org2", "Donor2_day70_dif1_org1",
    "Donor1_day70_dif2_org2", "Donor1_day70_dif2_org1", "Donor1_day70_dif1_org2",
    "Donor1_day70_dif1_org1", "Donor2_day40_dif2_org2", "Donor2_day40_dif2_org1",
    "Donor2_day40_dif1_org2", "Donor2_day40_dif1_org1", "Donor1_day40_dif2_org2",
    "Donor1_day40_dif2_org1", "Donor1_day40_dif1_org2", "Donor1_day40_dif1_org1",
]

rows = [
    [1],
    [1, 0.25],
    [1, 0.39, 0.15],
    [1, 0.59, 0.56, 0.26],
    [1, 0.46, 0.29, 0.6, 0.23],
    [1, 0.47, 0.58, 0.43, 0.47, 0.29],
    [1, 0.48, 0.36, 0.44, 0.32, 0.4, 0.52],
    [1, 0.25, 0.34, 0.23, 0.45, 0.81, 0.35, 0.1],
    [1, 0.26, 0.38, 0.55, 0.31, 0.6, 0.38, 0.34, 0.23],
    [1, 0.78, 0.25, 0.33, 0.47, 0.28, 0.58, 0.38, 0.31, 0.23],
    [1, 0.42, 0.44, 0.53, 0.37, 0.48, 0.37, 0.72, 0.6, 0.58, 0.19],
    [1, 0.63, 0.44, 0.47, 0.4, 0.51, 0.57, 0.57, 0.64, 0.46, 0.8, 0.38],
    [1, 0.31, 0.32, 0.34, 0.44, 0.16, 0.31, 0.41, 0.23, 0.44, 0.28, 0.24, 0.12],
    [1, 0.45, 0.5, 0.54, 0.49, 0.51, 0.29, 0.39, 0.48, 0.34, 0.69, 0.42, 0.42, 0.28],
    [1, 0.4, 0.33, 0.69, 0.64, 0.32, 0.37, 0.4, 0.37, 0.49, 0.47, 0.56, 0.43, 0.68, 0.17],
    [1, 0.49, 0.69, 0.42, 0.54, 0.7, 0.55, 0.57, 0.39, 0.35, 0.48, 0.34, 0.78, 0.52, 0.46, 0.2],
    [1, 0.18, 0.28, 0.18, 0.09, 0.37, 0.18, 0.16, 0.18, 0.07, 0.18, 0.33, 0.76, 0.29, 0.14, 0.42, 0.04],
    [1, 0.26, 0.55, 0.42, 0.47, 0.35, 0.53, 0.5, 0.74, 0.68, 0.42, 0.41, 0.59, 0.39, 0.63, 0.54, 0.4, 0.26],
    [1, 0.22, 0.68, 0.19, 0.33, 0.18, 0.09, 0.4, 0.18, 0.13, 0.16, 0.07, 0.14, 0.27, 0.62, 0.29, 0.13, 0.48, np.nan],
    [1, 0.15, 0.47, 0.16, 0.43, 0.47, 0.36, 0.25, 0.47, 0.53, 0.33, 0.36, 0.55, 0.35, 0.47, 0.31, 0.5, 0.58, 0.4, 0.17],
    [1, 0.4, 0.25, 0.5, 0.25, 0.7, 0.46, 0.85, 0.46, 0.58, 0.57, 0.5, 0.53, 0.32, 0.45, 0.53, 0.42, 0.74, 0.45, 0.5, 0.33],
    [1, 0.52, 0.4, 0.17, 0.78, 0.19, 0.56, 0.37, 0.5, 0.36, 0.49, 0.46, 0.85, 0.83, 0.32, 0.38, 0.55, 0.32, 0.61, 0.44, 0.35, 0.25],
    [1, 0.46, 0.6, 0.36, 0.16, 0.46, 0.2, 0.44, 0.39, 0.56, 0.49, 0.51, 0.41, 0.41, 0.48, 0.26, 0.58, 0.61, 0.38, 0.51, 0.35, 0.4, 0.51],
    [1, 0.38, 0.43, 0.53, 0.53, 0.17, 0.46, 0.16, 0.67, 0.65, 0.49, 0.3, 0.61, 0.88, 0.39, 0.41, 0.5, 0.34, 0.45, 0.34, 0.66, 0.54, 0.57, 0.17],
]

n = len(samples)
assert len(rows) == n
for i, r in enumerate(rows):
    assert len(r) == i + 1, f"row {i} ({samples[i]}) has {len(r)} values, expected {i+1}"

mat = pd.DataFrame(np.nan, index=samples, columns=samples)
for i, r in enumerate(rows):
    for k, v in enumerate(r):
        j = i - k  # sample[i] vs sample[i-k]
        mat.iloc[i, j] = v
        mat.iloc[j, i] = v

import sys

out = sys.argv[1] if len(sys.argv) > 1 else "Lab_Proprietary_Similarity_Matrix.csv"
mat.to_csv(out)
n_nan = mat.isna().sum().sum() - 0  # diagonal always filled here
print("Saved:", out)
print("Matrix shape:", mat.shape)
off = mat.values[~np.eye(n, dtype=bool)]
print("Off-diagonal NaN count:", np.isnan(off).sum(), "of", len(off))
