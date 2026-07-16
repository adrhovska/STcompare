#!/bin/bash

## User settings 
source "$(conda info --base)/etc/profile.d/conda.sh"
activate_env() {
  local env_name="$1"
  conda activate "$env_name"
  local status=$?
  if [[ "$status" -ne 0 ]]; then
    echo "Could not activate conda environment: $env_name" >&2
    exit 1
  fi
}
deactivate_env() {
  conda deactivate
}

PY_ENV="python_env"
R_ENV="r_env"
SKIP_LANDMARK_PICKING=false
MODE="single_pair"  
STCOMPARE_SCRIPT_NAME="STcompare.R"   # override with --stcompare_script to use a different R script (e.g. STcompareTissueClusters.r)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(pwd)"
SAMPLES_FILE="samples.txt"
REFERENCE_NAME=""

# help message
usage() {
  cat << EOF
Usage (single pair, original behaviour):
  $0 --source_dir <path> --reference_dir <path> \\
     --sample_aligned <name> --sample_reference <name> \\
     [--project_dir <path>] [--script_dir <path>] \\
     [--py_env <env>] [--r_env <env>] \\
     [--counts1 <path>] [--counts2 <path>] [--spatial2 <path>] \\
     [--skip_landmark_picking]

Usage (star approach, batch across all samples):
  $0 --mode star_landmarks  --samples_file samples.txt --reference_name <name> [--project_dir <path>]
  $0 --mode star_alignment  --samples_file samples.txt --reference_name <name> [--project_dir <path>]
  $0 --mode star_compare    --samples_file samples.txt --reference_name <name> [--project_dir <path>]
  $0 --mode star_all        --samples_file samples.txt --reference_name <name> [--project_dir <path>]
      (runs star_landmarks, star_alignment, star_compare in sequence)

Options:
  --mode                 single_pair (default) | star_landmarks | star_alignment | star_compare | star_all
  --samples_file         Tab-separated file: sample_name<TAB>space_ranger_outs_dir (for star modes)
  --reference_name        Name of the sample to use as the shared registration hub (for star modes)
  --source_dir            Source Space Ranger outs directory (single_pair mode)
  --reference_dir         Reference Space Ranger outs directory (single_pair mode)
  --sample_aligned        Name of sample being aligned (single_pair mode)
  --sample_reference      Name of reference sample (single_pair mode)
  --project_dir           Project/output directory
  --script_dir            Directory containing scripts
  --py_env                Conda Python environment
  --r_env                 Conda R environment
  --counts1               Source counts h5 file (single_pair mode)
  --counts2               Reference counts h5 file (single_pair mode)
  --spatial2              Reference spatial directory (single_pair mode)
  --skip_landmark_picking Skip LandmarkPicker.py, reuse existing landmark CSVs (single_pair mode)
  --stcompare_script       Name of the R script (in --script_dir) to run for the STcompare step (default: STcompare.R)
  --help                  Show help message
EOF
}

need_value() {
  if [[ $# -lt 2 || "$2" == -* ]]; then
    echo "Missing value for $1" >&2
    usage
    exit 1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) need_value "$@"; MODE="$2"; shift 2 ;;
    --samples_file) need_value "$@"; SAMPLES_FILE="$2"; shift 2 ;;
    --reference_name) need_value "$@"; REFERENCE_NAME="$2"; shift 2 ;;
    --source_dir) need_value "$@"; SOURCE_DIR="$2"; shift 2 ;;
    --reference_dir) need_value "$@"; REFERENCE_DIR="$2"; shift 2 ;;
    --sample_aligned) need_value "$@"; SAMPLE_ALIGNED="$2"; shift 2 ;;
    --sample_reference) need_value "$@"; SAMPLE_REFERENCE="$2"; shift 2 ;;
    --project_dir) need_value "$@"; PROJECT_DIR="$2"; shift 2 ;;
    --script_dir) need_value "$@"; SCRIPT_DIR="$2"; shift 2 ;;
    --py_env) need_value "$@"; PY_ENV="$2"; shift 2 ;;
    --r_env) need_value "$@"; R_ENV="$2"; shift 2 ;;
    --counts1) need_value "$@"; COUNTS1="$2"; shift 2 ;;
    --counts2) need_value "$@"; COUNTS2="$2"; shift 2 ;;
    --spatial2) need_value "$@"; SPATIAL2="$2"; shift 2 ;;
    --skip_landmark_picking) SKIP_LANDMARK_PICKING=true; shift 1 ;;
    --stcompare_script) need_value "$@"; STCOMPARE_SCRIPT_NAME="$2"; shift 2 ;;
    --help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

LANDMARK_PICKER="${SCRIPT_DIR}/scripts/landmarking/LandmarkPicker.py"
STALIGN_SCRIPT="${SCRIPT_DIR}/scripts/alignment/STalignCode.py"
QC_PLOTS_SCRIPT="${SCRIPT_DIR}/scripts/alignment/STalignQC.py"
STCOMPARE_SCRIPT="${SCRIPT_DIR}/scripts/comparison/${STCOMPARE_SCRIPT_NAME}"

need_file() { if [[ ! -f "$1" ]]; then echo "Missing file: $1"; exit 1; fi }
need_dir()  { if [[ ! -d "$1" ]]; then echo "Missing directory: $1"; exit 1; fi }
require_arg() {
  local value="$1" name="$2"
  if [[ -z "$value" ]]; then
    echo "Missing required argument: $name" >&2
    usage
    exit 1
  fi
}

## 1. Core single-pair workflow (landmarks -> STalign -> STcompare)
--
run_single_pair() {
  local src_dir="$1" ref_dir="$2" sample_aligned="$3" sample_reference="$4"
  local project_dir="$5" skip_landmarks="$6"
  local counts1="${7:-}" counts2="${8:-}" spatial2="${9:-}"

  local source_image="${src_dir}/spatial/tissue_hires_image.png"
  local reference_image="${ref_dir}/spatial/tissue_hires_image.png"
  local source_pos="${src_dir}/spatial/tissue_positions.csv"
  local reference_pos="${ref_dir}/spatial/tissue_positions.csv"
  local source_scale="${src_dir}/spatial/scalefactors_json.json"
  local reference_scale="${ref_dir}/spatial/scalefactors_json.json"

  counts1="${counts1:-${src_dir}/filtered_feature_bc_matrix.h5}"
  counts2="${counts2:-${ref_dir}/filtered_feature_bc_matrix.h5}"
  spatial2="${spatial2:-${ref_dir}/spatial}"

  local align_pair_name="${sample_aligned}_aligned_to_${sample_reference}"
  local landmark_pair_name="${sample_aligned}_paired_to_${sample_reference}"
  local run_dir="${project_dir}/STcompare_outputs/${align_pair_name}"
  local landmark_dir="${run_dir}/landmarks/${landmark_pair_name}"
  local stalign_outdir="${run_dir}/STalign"
  local stcompare_outdir="${run_dir}/STcompare"
  local points1="${landmark_dir}/${sample_aligned}_points.csv"
  local points2="${landmark_dir}/${sample_reference}_points.csv"
  local aligned_pos="${stalign_outdir}/${sample_aligned}_aligned_to_${sample_reference}_barcodes.csv"

  need_dir "$src_dir"; need_dir "$ref_dir"; need_dir "$spatial2"
  mkdir -p "$run_dir" "$landmark_dir" "$stalign_outdir" "$stcompare_outdir"

  if [[ "$skip_landmarks" == true ]]; then
    echo "Step 1: Skipping LandmarkPicker.py, reusing existing landmarks for $align_pair_name"
    need_file "$points1"; need_file "$points2"
  else
    echo "Step 1: Running LandmarkPicker.py for $align_pair_name"
    activate_env "$PY_ENV"
    python "$LANDMARK_PICKER" \
      --image1 "$source_image" --image2 "$reference_image" \
      --sample_aligned "$sample_aligned" --sample_reference "$sample_reference" \
      --project_dir "$run_dir"
    deactivate_env
    need_file "$points1"; need_file "$points2"
  fi

  echo "Step 2: Running STalignCode.py for $align_pair_name"
  activate_env "$PY_ENV"
  PYTHONPATH="$SCRIPT_DIR/scripts/alignment:${PYTHONPATH:-}" python "$STALIGN_SCRIPT" \
    --image1 "$source_image" --image2 "$reference_image" \
    --pos1 "$source_pos" --pos2 "$reference_pos" \
    --scale1 "$source_scale" --scale2 "$reference_scale" \
    --points1 "$points1" --points2 "$points2" \
    --sample_aligned "$sample_aligned" --sample_reference "$sample_reference" \
    --project_dir "$run_dir" --outdir "$stalign_outdir" \
    --alignment_method stalign
  deactivate_env
  need_file "$aligned_pos"

  echo "Step 3: Running STcompare.R for $align_pair_name"
  activate_env "$R_ENV"
  Rscript "$STCOMPARE_SCRIPT" \
    --counts1 "$counts1" --counts2 "$counts2" \
    --pos1 "$aligned_pos" --type1 aligned \
    --pos2 "$spatial2" --type2 visium \
    --outdir "$stcompare_outdir" \
    --sample_aligned "$sample_aligned" --sample_reference "$sample_reference"
  deactivate_env

  echo "Pair complete: $align_pair_name"
}
# load samples.txt into NAMES[]/DIRS[] arrays for star modes

load_samples() {
  need_file "$SAMPLES_FILE"
  NAMES=(); DIRS=()
  while IFS=$'\t' read -r name dir; do
    NAMES+=("$name"); DIRS+=("$dir")
  done < "$SAMPLES_FILE"
}

get_dir_for_name() {
  local target="$1"
  for i in "${!NAMES[@]}"; do
    [[ "${NAMES[$i]}" == "$target" ]] && { echo "${DIRS[$i]}"; return; }
  done
}

# star mode: landmarks (each sample vs the reference, once each)

run_star_landmarks() {
  require_arg "$REFERENCE_NAME" "--reference_name"
  load_samples
  local ref_dir; ref_dir="$(get_dir_for_name "$REFERENCE_NAME")"
  require_arg "$ref_dir" "reference sample directory (not found in $SAMPLES_FILE)"

  for i in "${!NAMES[@]}"; do
    local name="${NAMES[$i]}" dir="${DIRS[$i]}"
    [[ "$name" == "$REFERENCE_NAME" ]] && continue

    local align_pair_name="${name}_aligned_to_${REFERENCE_NAME}"
    local landmark_pair_name="${name}_paired_to_${REFERENCE_NAME}"
    local run_dir="${PROJECT_DIR}/STcompare_outputs/${align_pair_name}"
    local points1="${run_dir}/landmarks/${landmark_pair_name}/${name}_points.csv"
    local points2="${run_dir}/landmarks/${landmark_pair_name}/${REFERENCE_NAME}_points.csv"

    if [[ -f "$points1" && -f "$points2" ]]; then
      echo "Skipping already-picked sample: $name"
      continue
    fi

    mkdir -p "$run_dir"
    echo ""
    "Pick landmarks: $name  vs  reference ($REFERENCE_NAME)"
    echo ""

    activate_env "$PY_ENV"
    python "$LANDMARK_PICKER" \
      --image1 "${dir}/spatial/tissue_hires_image.png" \
      --image2 "${ref_dir}/spatial/tissue_hires_image.png" \
      --sample_aligned "$name" --sample_reference "$REFERENCE_NAME" \
      --project_dir "$run_dir"
    deactivate_env
  done
}

# star mode: alignment (register every sample into reference frame)
run_star_alignment() {
  require_arg "$REFERENCE_NAME" "--reference_name"
  load_samples
  local ref_dir; ref_dir="$(get_dir_for_name "$REFERENCE_NAME")"
  require_arg "$ref_dir" "reference sample directory (not found in $SAMPLES_FILE)"

  local pair_log="${PROJECT_DIR}/star_alignment_completed.txt"
  mkdir -p "$PROJECT_DIR"; touch "$pair_log"

  for i in "${!NAMES[@]}"; do
    local name="${NAMES[$i]}" dir="${DIRS[$i]}"
    [[ "$name" == "$REFERENCE_NAME" ]] && continue
    if grep -qx "$name" "$pair_log"; then
      echo "Skipping already-aligned sample: $name"
      continue
    fi

    local align_pair_name="${name}_aligned_to_${REFERENCE_NAME}"
    local landmark_pair_name="${name}_paired_to_${REFERENCE_NAME}"
    local run_dir="${PROJECT_DIR}/STcompare_outputs/${align_pair_name}"
    local points1="${run_dir}/landmarks/${landmark_pair_name}/${name}_points.csv"
    local points2="${run_dir}/landmarks/${landmark_pair_name}/${REFERENCE_NAME}_points.csv"
    local stalign_outdir="${run_dir}/STalign"

    if [[ ! -f "$points1" || ! -f "$points2" ]]; then
      echo "WARNING: no landmarks for $name, skipping (run --mode star_landmarks first)"
      continue
    fi

    mkdir -p "$stalign_outdir"
    echo "Aligning $name into reference frame ($REFERENCE_NAME)"
    
    activate_env "$PY_ENV"
    PYTHONPATH="$SCRIPT_DIR/scripts/alignment:${PYTHONPATH:-}" python "$STALIGN_SCRIPT" \
      --image1 "${dir}/spatial/tissue_hires_image.png" \
      --image2 "${ref_dir}/spatial/tissue_hires_image.png" \
      --pos1 "${dir}/spatial/tissue_positions.csv" \
      --pos2 "${ref_dir}/spatial/tissue_positions.csv" \
      --scale1 "${dir}/spatial/scalefactors_json.json" \
      --scale2 "${ref_dir}/spatial/scalefactors_json.json" \
      --points1 "$points1" --points2 "$points2" \
      --sample_aligned "$name" --sample_reference "$REFERENCE_NAME" \
      --project_dir "$run_dir" --outdir "$stalign_outdir" \
      --alignment_method stalign
    deactivate_env

    echo "$name" >> "$pair_log"
  done
}

# star mode: compare (every unique pair, using shared reference frame)

get_sample_info() {
  local sample_name="$1" dir
  dir="$(get_dir_for_name "$sample_name")"
  if [[ "$sample_name" == "$REFERENCE_NAME" ]]; then
    echo "${dir}/filtered_feature_bc_matrix.h5|${dir}/spatial|visium"
  else
    local align_pair="${sample_name}_aligned_to_${REFERENCE_NAME}"
    local aligned_csv="${PROJECT_DIR}/STcompare_outputs/${align_pair}/STalign/${sample_name}_aligned_to_${REFERENCE_NAME}_barcodes.csv"
    echo "${dir}/filtered_feature_bc_matrix.h5|${aligned_csv}|aligned"
  fi
}

run_star_compare() {
  require_arg "$REFERENCE_NAME" "--reference_name"
  load_samples
  local n=${#NAMES[@]}
  local pair_log="${PROJECT_DIR}/star_compare_completed.txt"
  mkdir -p "$PROJECT_DIR"; touch "$pair_log"

  activate_env "$R_ENV"
  for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do
      local s1="${NAMES[$i]}" s2="${NAMES[$j]}"
      local pair_tag="${s1}_vs_${s2}"

      if grep -qx "$pair_tag" "$pair_log"; then
        echo "Skipping already-compared pair: $pair_tag"
        continue
      fi

      IFS='|' read -r counts1 pos1 type1 <<< "$(get_sample_info "$s1")"
      IFS='|' read -r counts2 pos2 type2 <<< "$(get_sample_info "$s2")"

      if [[ "$type1" == "aligned" && ! -f "$pos1" ]]; then
        echo "WARNING: missing aligned coords for $s1, skipping $pair_tag"
        continue
      fi
      if [[ "$type2" == "aligned" && ! -f "$pos2" ]]; then
        echo "WARNING: missing aligned coords for $s2, skipping $pair_tag"
        continue
      fi

      local stcompare_outdir="${PROJECT_DIR}/STcompare_pairwise/${pair_tag}"
      mkdir -p "$stcompare_outdir"

      echo "Comparing: $pair_tag"

      Rscript "$STCOMPARE_SCRIPT" \
        --counts1 "$counts1" --counts2 "$counts2" \
        --pos1 "$pos1" --type1 "$type1" \
        --pos2 "$pos2" --type2 "$type2" \
        --outdir "$stcompare_outdir" \
        --sample_aligned "$s1" --sample_reference "$s2"

      echo "$pair_tag" >> "$pair_log"
    done
  done
  deactivate_env
}

case "$MODE" in
  single_pair)
    require_arg "${SOURCE_DIR:-}" "--source_dir"
    require_arg "${REFERENCE_DIR:-}" "--reference_dir"
    require_arg "${SAMPLE_ALIGNED:-}" "--sample_aligned"
    require_arg "${SAMPLE_REFERENCE:-}" "--sample_reference"
    run_single_pair "$SOURCE_DIR" "$REFERENCE_DIR" "$SAMPLE_ALIGNED" "$SAMPLE_REFERENCE" \
      "$PROJECT_DIR" "$SKIP_LANDMARK_PICKING" "${COUNTS1:-}" "${COUNTS2:-}" "${SPATIAL2:-}"
    ;;
  star_landmarks) run_star_landmarks ;;
  star_alignment) run_star_alignment ;;
  star_compare)   run_star_compare ;;
  star_all)
    run_star_landmarks
    run_star_alignment
    run_star_compare
    ;;
  *)
    echo "Unknown --mode: $MODE" >&2
    usage
    exit 1
    ;;
esac

echo "Complete :D"