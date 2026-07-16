#!/usr/bin/env bash

# Runs STcompare_with_tissue_clusters.R for every unique pair among your 24
# organoid samples, then builds the per-tissue-type-cluster heatmaps

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=true
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="/path/to/your/organoid/data"
ALIGNED_ROOT="/path/to/your/aligned/positions"
OUT_ROOT="/path/to/STcompare_all_runs"
SAMPLES=(
  Donor1_day40_dif1_org1  Donor1_day40_dif1_org2
  Donor1_day40_dif2_org1  Donor1_day40_dif2_org2
  Donor2_day40_dif1_org1  Donor2_day40_dif1_org2
  Donor2_day40_dif2_org1  Donor2_day40_dif2_org2
  Donor1_day70_dif1_org1  Donor1_day70_dif1_org2
  Donor1_day70_dif2_org1  Donor1_day70_dif2_org2
  Donor2_day70_dif1_org1  Donor2_day70_dif1_org2
  Donor2_day70_dif2_org1  Donor2_day70_dif2_org2
  Donor1_day120_dif1_org1 Donor1_day120_dif1_org2
  Donor1_day120_dif2_org1 Donor1_day120_dif2_org2
  Donor2_day120_dif1_org1 Donor2_day120_dif1_org2
  Donor2_day120_dif2_org1 Donor2_day120_dif2_org2
)

RES=20         
THREADS=4      
SCALE=hires    

counts_path()  { echo "$DATA_ROOT/$1/filtered_feature_bc_matrix.h5"; }
visium_path()  { echo "$DATA_ROOT/$1/spatial"; }
aligned_path() { echo "$ALIGNED_ROOT/${1}_to_${2}_aligned.csv"; }

mkdir -p "$OUT_ROOT"

n=${#SAMPLES[@]}
total_pairs=$(( n * (n - 1) / 2 ))
echo "Planning $total_pairs pairs across $n samples"

count=0
skipped=0
failed=0

for (( i=0; i<n; i++ )); do
  for (( j=i+1; j<n; j++ )); do
    sampleA="${SAMPLES[$i]}"
    sampleB="${SAMPLES[$j]}"
    count=$((count + 1))

    counts1="$(counts_path "$sampleA")"
    counts2="$(counts_path "$sampleB")"
    pos1="$(aligned_path "$sampleA" "$sampleB")"
    pos2="$(visium_path "$sampleB")"
    outdir="$OUT_ROOT/${sampleA}_vs_${sampleB}"

    cmd=(Rscript "$SCRIPT_DIR/STcompare_with_tissue_clusters.R"
      --counts1 "$counts1"
      --counts2 "$counts2"
      --pos1 "$pos1" --type1 aligned
      --pos2 "$pos2" --type2 visium
      --outdir "$outdir"
      --scale "$SCALE" --res "$RES" --threads "$THREADS"
      --sample_aligned "$sampleA" --sample_reference "$sampleB"
    )

    echo "[$count/$total_pairs] $sampleA vs $sampleB"

    if $DRY_RUN; then
      printf '  %q' "${cmd[@]}"; echo
      continue
    fi

    # skip pairs that are missing required input files instead of failing outright
    if [[ ! -f "$counts1" || ! -f "$counts2" || ! -f "$pos1" ]]; then
      echo "  SKIP: missing input file(s) for this pair"
      echo "    counts1: $counts1 $( [[ -f "$counts1" ]] && echo OK || echo MISSING )"
      echo "    counts2: $counts2 $( [[ -f "$counts2" ]] && echo OK || echo MISSING )"
      echo "    pos1(aligned): $pos1 $( [[ -f "$pos1" ]] && echo OK || echo MISSING )"
      skipped=$((skipped + 1))
      continue
    fi

    if "${cmd[@]}"; then
      echo "  done -> $outdir"
    else
      echo "  FAILED: $sampleA vs $sampleB"
      failed=$((failed + 1))
    fi
  done
done

echo ""
echo "Ran $((count - skipped)) pairs, skipped $skipped (missing files), $failed failed."

if $DRY_RUN; then
  echo "Dry run only - no heatmaps built."
  exit 0
fi

echo "Building tissue-type cluster heatmaps from all completed pairs..."
Rscript "$SCRIPT_DIR/build_tissue_cluster_heatmaps.R" \
  --stcompare_root "$OUT_ROOT" \
  --outdir "$OUT_ROOT/Cluster_Heatmaps"

echo "Done. See $OUT_ROOT/Cluster_Heatmaps for the per-cluster similarity heatmaps."