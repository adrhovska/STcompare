"""
Builds, from the STcompare pairwise output folders:
  1. the overall similarity matrix/heatmap (mean_percent_similarity), reordered
     to match the reference heatmap ordering used in SimilarityHeatmap.r
  2. three cluster-level correlation matrices/heatmaps (progenitor_genes,
     maturation_genes, patterning_genes), using correlationCoef from
     Cluster_Level_Results.csv

Sample ordering matches SimilarityHeatmap.r's order_samples(): sorted by
day (120,70,40), then donor (Donor2,Donor1), then dif (dif2,dif1), then
org (org2,org1) -- i.e. descending on every component. This matches the
lab's reference "all_organoids_spatial_similarity_heatmap.pdf" ordering.
"""
import glob
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "."
os.makedirs(OUT_DIR, exist_ok=True)

DAY_ORDER = ["day120", "day70", "day40"]
DONOR_ORDER = ["Donor2", "Donor1"]
DIF_ORDER = ["dif2", "dif1"]
ORG_ORDER = ["org2", "org1"]
PAT = re.compile(r"^(Donor\d+)_(day\d+)_(dif\d+)_(org\d+)$")


def order_samples(names):
    def key(s):
        m = PAT.match(s)
        if not m:
            raise ValueError(f"Unexpected sample name: {s}")
        donor, day, dif, org = m.groups()
        return (DAY_ORDER.index(day), DONOR_ORDER.index(donor),
                DIF_ORDER.index(dif), ORG_ORDER.index(org))
    return sorted(names, key=key)


pair_dirs = sorted(glob.glob(os.path.join(ROOT, "*_vs_*")))
pair_dirs = [d for d in pair_dirs if os.path.isdir(d)]

overall_rows, cluster_rows = [], []
missing = []
for d in pair_dirs:
    base = os.path.basename(d)
    a, b = base.split("_vs_")
    if a == b:
        continue  # known self-comparison artefact
    ov = os.path.join(d, "Results", "Overall_Similarity.csv")
    cl = os.path.join(d, "Results", "Cluster_Level_Results.csv")
    if os.path.isfile(ov):
        overall_rows.append(pd.read_csv(ov).iloc[0])
    else:
        missing.append(ov)
    if os.path.isfile(cl):
        cdf = pd.read_csv(cl, index_col=0)
        cdf["sample_aligned"] = a
        cdf["sample_reference"] = b
        cluster_rows.append(cdf.reset_index().rename(columns={"index": "cluster"}))
    else:
        missing.append(cl)

if missing:
    print(f"{len(missing)} missing files (excluding known self-comparison), e.g.:")
    for m in missing[:5]:
        print("  ", m)

overall = pd.DataFrame(overall_rows).reset_index(drop=True)
overall["pair_key"] = overall.apply(lambda r: tuple(sorted([r.sample_aligned, r.sample_reference])), axis=1)
overall = overall.drop_duplicates(subset="pair_key", keep="first").reset_index(drop=True)

cluster_long = pd.concat(cluster_rows, ignore_index=True)
cluster_long["pair_key"] = cluster_long.apply(lambda r: tuple(sorted([r.sample_aligned, r.sample_reference])), axis=1)
cluster_long = cluster_long.drop_duplicates(subset=["pair_key", "cluster"], keep="first").reset_index(drop=True)

samples = order_samples(set(overall.sample_aligned) | set(overall.sample_reference))
n = len(samples)
print(f"Samples: {n}, expected pairs: {n*(n-1)//2}, overall rows: {len(overall)}")


def build_matrix(df, value_col, samples, diag_value=1.0):
    mat = pd.DataFrame(np.nan, index=samples, columns=samples)
    for _, r in df.iterrows():
        a, b, v = r.sample_aligned, r.sample_reference, r[value_col]
        mat.loc[a, b] = v
        mat.loc[b, a] = v
    np.fill_diagonal(mat.values, diag_value)
    return mat


def plot_seq_heatmap(mat, title, out_pdf, label="Similarity"):
    off = mat.values[~np.eye(len(mat), dtype=bool)]
    off = off[~np.isnan(off)]
    vmax = off.max()
    plt.figure(figsize=(11, 10))
    sns.heatmap(mat, vmin=0, vmax=vmax, cmap="magma", square=True,
                annot=True, fmt=".2f", annot_kws={"size": 6},
                cbar_kws={"label": label}, linewidths=0.4, linecolor="white")
    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def plot_div_heatmap(mat, title, out_pdf, label="Correlation coefficient"):
    plt.figure(figsize=(11, 10))
    sns.heatmap(mat, vmin=-1, vmax=1, cmap="RdBu_r", center=0, square=True,
                annot=True, fmt=".2f", annot_kws={"size": 6},
                cbar_kws={"label": label}, linewidths=0.4, linecolor="white")
    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


# 1. overall similarity matrix, reordered
overall_mat = build_matrix(overall, "mean_percent_similarity", samples, diag_value=1.0)
overall_mat.to_csv(os.path.join(OUT_DIR, "All_Pairwise_Similarity_Matrix.csv"))
overall.drop(columns="pair_key").to_csv(os.path.join(OUT_DIR, "All_Pairwise_Similarity_Long.csv"), index=False)
plot_seq_heatmap(overall_mat, "Overall spatial similarity (all genes)",
                  os.path.join(OUT_DIR, "spatial_similarity_heatmap.pdf"))

# 2. cluster-level correlation matrices
cluster_names = ["progenitor_genes", "maturation_genes", "patterning_genes"]
for cname in cluster_names:
    sub = cluster_long[cluster_long["cluster"] == cname]
    n_valid = sub["correlationCoef"].notna().sum()
    print(f"{cname}: {len(sub)} pairs, {n_valid} with a valid correlationCoef")
    mat = build_matrix(sub, "correlationCoef", samples, diag_value=1.0)
    mat.to_csv(os.path.join(OUT_DIR, f"Cluster_{cname}_Correlation_Matrix.csv"))
    plot_div_heatmap(mat, f"Cluster-level spatial correlation: {cname}",
                      os.path.join(OUT_DIR, f"Cluster_{cname}_correlation_heatmap.pdf"))

cluster_long.drop(columns="pair_key").to_csv(
    os.path.join(OUT_DIR, "All_Cluster_Level_Results_Long.csv"), index=False
)

print(f"\nDone. Outputs written to {OUT_DIR}")
