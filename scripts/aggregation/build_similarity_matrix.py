"""
Aggregate STcompare pairwise outputs into a single sample x sample
similarity matrix + heatmap.

Reads: <PAIR_DIR>/<sample_aligned>_vs_<sample_reference>/Results/Overall_Similarity.csv
Writes: similarity_matrix.csv, similarity_matrix_long.csv, OwnHeatmap.pdf

Handles two known data-quality issues found in this dataset:
  1. A self-comparison folder (Donor1_day40_dif2_org2_vs_Donor1_day40_dif2_org2)
     - excluded, since a sample cannot be meaningfully compared to itself.
  2. 10 sample pairs that were run in both directions (A_vs_B and B_vs_A),
     all involving Donor1_day40_dif2_org2 - both directions give identical
     mean_percent_similarity values (confirmed), so the duplicate is dropped
     and only one direction is kept.
After cleanup: 24 samples, 276 unique pairs (= C(24,2)), matching the
pairwise design described in the thesis.
"""
import glob
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "."

os.makedirs(OUT_DIR, exist_ok=True)

rows = []
missing = []
pair_dirs = sorted(glob.glob(os.path.join(ROOT, "*_vs_*")))

for d in pair_dirs:
    csv_path = os.path.join(d, "Results", "Overall_Similarity.csv")
    if not os.path.isfile(csv_path):
        missing.append(d)
        continue
    df = pd.read_csv(csv_path)
    rows.append(df.iloc[0])

raw = pd.DataFrame(rows).reset_index(drop=True)
print(f"Pair folders found: {len(pair_dirs)}")
print(f"Missing Overall_Similarity.csv: {len(missing)}")
for m in missing:
    print("  MISSING:", m)
print(f"Rows parsed: {len(raw)}")

# canonical (unordered) pair key
raw["pair_key"] = raw.apply(
    lambda r: tuple(sorted([r["sample_aligned"], r["sample_reference"]])), axis=1
)

# 1. drop self-comparisons
self_comp = raw[raw["sample_aligned"] == raw["sample_reference"]]
if len(self_comp):
    print("\nExcluding self-comparison rows:")
    print(self_comp[["sample_aligned", "sample_reference", "mean_percent_similarity"]])
raw = raw[raw["sample_aligned"] != raw["sample_reference"]].copy()

# 2. check direction-duplicates agree, then dedup
dup_keys = raw["pair_key"][raw["pair_key"].duplicated(keep=False)].unique()
if len(dup_keys):
    print(f"\n{len(dup_keys)} pairs run in both directions - checking agreement:")
    for k in dup_keys:
        sub = raw[raw["pair_key"] == k]
        vals = sub["mean_percent_similarity"].round(10).unique()
        status = "OK (identical)" if len(vals) == 1 else "MISMATCH!!"
        print(f"  {k}: {list(sub['mean_percent_similarity'])} -> {status}")

clean = raw.drop_duplicates(subset="pair_key", keep="first").reset_index(drop=True)
print(f"\nRows after removing self-comparison + direction duplicates: {len(clean)}")

samples = sorted(set(clean["sample_aligned"]) | set(clean["sample_reference"]))
n = len(samples)
print(f"Unique samples: {n}  (expected pairs = {n*(n-1)//2})")

mat = pd.DataFrame(np.nan, index=samples, columns=samples)
for _, r in clean.iterrows():
    a, b, v = r["sample_aligned"], r["sample_reference"], r["mean_percent_similarity"]
    mat.loc[a, b] = v
    mat.loc[b, a] = v
np.fill_diagonal(mat.values, 1.0)  # hardcoded self-comparison, as in existing figure caption

off_diag = mat.values[~np.eye(n, dtype=bool)]
off_diag = off_diag[~np.isnan(off_diag)]
print(f"\nOff-diagonal stats: n={len(off_diag)}  min={off_diag.min():.3f}  "
      f"max={off_diag.max():.3f}  mean={off_diag.mean():.3f}  sd={off_diag.std(ddof=1):.3f}")

nan_cells = np.isnan(mat.values) & ~np.eye(n, dtype=bool)
if nan_cells.any():
    print(f"\nWARNING: {nan_cells.sum()//2} sample pairs missing from matrix:")
    idx = np.where(np.triu(nan_cells, k=1))
    for i, j in zip(*idx):
        print(f"  {samples[i]}  vs  {samples[j]}")

mat.to_csv(os.path.join(OUT_DIR, "similarity_matrix.csv"))
clean[["sample_aligned", "sample_reference", "mean_percent_similarity", "n_genes_evaluated"]].to_csv(
    os.path.join(OUT_DIR, "similarity_matrix_long.csv"), index=False
)

# heatmap, colour capped at max off-diagonal value (as in existing figure caption)
vmax = off_diag.max()
plt.figure(figsize=(10, 9))
sns.heatmap(
    mat, vmin=0, vmax=vmax, cmap="viridis",
    square=True, cbar_kws={"label": "Mean spatial similarity"},
    xticklabels=True, yticklabels=True,
)
plt.title("Pairwise spatial similarity (STcompare)")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "OwnHeatmap.pdf"))
plt.savefig(os.path.join(OUT_DIR, "OwnHeatmap.png"), dpi=200)
print(f"\nSaved: similarity_matrix.csv, similarity_matrix_long.csv, OwnHeatmap.pdf/png -> {OUT_DIR}")
