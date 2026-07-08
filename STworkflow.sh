#!/bin/bash

## User settings 
# conda environments setup and activation
# setup environment names based on your local conda environment names
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

# project folder 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(pwd)"

# help message
usage() {
  cat << EOF
Usage: $0 \\
  --source_dir <path> \\
  --reference_dir <path> \\
  --sample_aligned <name> \\
  --sample_reference <name> \\
  [--project_dir <path>] \\
  [--script_dir <path>] \\
  [--py_env <env>] \\
  [--r_env <env>] \\
  [--counts1 <path>] \\
  [--counts2 <path>] \\
  [--spatial2 <path>]

Required:
  --source_dir          Source Space Ranger outs directory
  --reference_dir       Reference Space Ranger outs directory
  --sample_aligned      Name of sample being aligned
  --sample_reference    Name of reference sample

Optional:
  --project_dir         Project/output directory
  --script_dir          Directory containing scripts
  --py_env              Conda Python environment
  --r_env               Conda R environment
  --counts1             Source counts h5 file
  --counts2             Reference counts h5 file
  --spatial2            Reference spatial directory
  --help                Show help message
EOF
}

need_value() {
  if [[ $# -lt 2 || "$2" == -* ]]; then
    echo "Missing value for $1" >&2
    usage
    exit 1
  fi
}

# parse arguments 
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source_dir)
      need_value "$@"
      SOURCE_DIR="$2"
      shift 2
      ;;
    --reference_dir)
      need_value "$@"
      REFERENCE_DIR="$2"
      shift 2
      ;;
    --sample_aligned)
      need_value "$@"
      SAMPLE_ALIGNED="$2"
      shift 2
      ;;
    --sample_reference)
      need_value "$@"
      SAMPLE_REFERENCE="$2"
      shift 2
      ;;
    --project_dir)
      need_value "$@"
      PROJECT_DIR="$2"
      shift 2
      ;;
    --script_dir)
      need_value "$@"
      SCRIPT_DIR="$2"
      shift 2
      ;;
    --py_env)
      need_value "$@"
      PY_ENV="$2"
      shift 2
      ;;
    --r_env)
      need_value "$@"
      R_ENV="$2"
      shift 2
      ;;
    --counts1)
      need_value "$@"
      COUNTS1="$2"
      shift 2
      ;;
    --counts2)
      need_value "$@"
      COUNTS2="$2"
      shift 2
      ;;
    --spatial2)
      need_value "$@"
      SPATIAL2="$2"
      shift 2
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# check required arguments
require_arg() {
  local value="$1"
  local name="$2"

  if [[ -z "$value" ]]; then
    echo "Missing required argument: $name" >&2
    usage
    exit 1
  fi
}
require_arg "${SOURCE_DIR:-}" "--source_dir"
require_arg "${REFERENCE_DIR:-}" "--reference_dir"
require_arg "${SAMPLE_ALIGNED:-}" "--sample_aligned"
require_arg "${SAMPLE_REFERENCE:-}" "--sample_reference"

# derived paths
LANDMARK_PICKER="${SCRIPT_DIR}/LandmarkPicker.py"
STALIGN_SCRIPT="${SCRIPT_DIR}/STalignCode.py"
STCOMPARE_SCRIPT="${SCRIPT_DIR}/STcompare.R"

SOURCE_IMAGE="${SOURCE_DIR}/spatial/tissue_hires_image.png"
REFERENCE_IMAGE="${REFERENCE_DIR}/spatial/tissue_hires_image.png"

SOURCE_POS="${SOURCE_DIR}/spatial/tissue_positions.csv"
REFERENCE_POS="${REFERENCE_DIR}/spatial/tissue_positions.csv"

SOURCE_SCALE="${SOURCE_DIR}/spatial/scalefactors_json.json"
REFERENCE_SCALE="${REFERENCE_DIR}/spatial/scalefactors_json.json"

COUNTS1="${COUNTS1:-${SOURCE_DIR}/filtered_feature_bc_matrix.h5}"
COUNTS2="${COUNTS2:-${REFERENCE_DIR}/filtered_feature_bc_matrix.h5}"
SPATIAL2="${SPATIAL2:-${REFERENCE_DIR}/spatial}"

ALIGN_PAIR_NAME="${SAMPLE_ALIGNED}_aligned_to_${SAMPLE_REFERENCE}"
LANDMARK_PAIR_NAME="${SAMPLE_ALIGNED}_paired_to_${SAMPLE_REFERENCE}"

# one shared output folder for the full workflow
RUN_DIR="${PROJECT_DIR}/STcompare_outputs/${ALIGN_PAIR_NAME}"

# outputs from each step
LANDMARK_DIR="${RUN_DIR}/landmarks/${LANDMARK_PAIR_NAME}"
STALIGN_OUTDIR="${RUN_DIR}/STalign"
STCOMPARE_OUTDIR="${RUN_DIR}/STcompare"

POINTS1="${LANDMARK_DIR}/${SAMPLE_ALIGNED}_points.csv"
POINTS2="${LANDMARK_DIR}/${SAMPLE_REFERENCE}_points.csv"

ALIGNED_POS="${STALIGN_OUTDIR}/${SAMPLE_ALIGNED}_aligned_to_${SAMPLE_REFERENCE}_barcodes.csv"

# helper checks for required files and directories
need_file() {
  if [[ ! -f "$1" ]]; then
    echo "Missing file: $1"
    exit 1
  fi
}
need_dir() {
  if [[ ! -d "$1" ]]; then
    echo "Missing directory: $1"
    exit 1
  fi
}

need_dir "$SOURCE_DIR"
need_dir "$REFERENCE_DIR"
need_dir "$SPATIAL2"

# create main run directory
mkdir -p "$RUN_DIR" "$LANDMARK_DIR" "$STALIGN_OUTDIR" "$STCOMPARE_OUTDIR"

## Workflow
# 1: LandmarkPicker.py
echo "Step 1: Running LandmarkPicker.py"
echo "Click matching landmarks"
activate_env "$PY_ENV"
python "$LANDMARK_PICKER" \
  --image1 "$SOURCE_IMAGE" \
  --image2 "$REFERENCE_IMAGE" \
  --sample_aligned "$SAMPLE_ALIGNED" \
  --sample_reference "$SAMPLE_REFERENCE" \
  --project_dir "$RUN_DIR"
deactivate_env

# check landmark files
need_file "$POINTS1"
need_file "$POINTS2"

# 2: STalignCode.py
echo "Step 2: Running STalignCode.py"
activate_env "$PY_ENV"
python "$STALIGN_SCRIPT" \
  --image1 "$SOURCE_IMAGE" \
  --image2 "$REFERENCE_IMAGE" \
  --pos1 "$SOURCE_POS" \
  --pos2 "$REFERENCE_POS" \
  --scale1 "$SOURCE_SCALE" \
  --scale2 "$REFERENCE_SCALE" \
  --points1 "$POINTS1" \
  --points2 "$POINTS2" \
  --sample_aligned "$SAMPLE_ALIGNED" \
  --sample_reference "$SAMPLE_REFERENCE" \
  --project_dir "$RUN_DIR" \
  --outdir "$STALIGN_OUTDIR" \
  --alignment_method stalign
deactivate_env

# check aligned coordinate file
need_file "$ALIGNED_POS"

# 3: STcompare.R
echo "Step 3: Running STcompare.R"
activate_env "$R_ENV"
Rscript "$STCOMPARE_SCRIPT" \
  --counts1 "$COUNTS1" \
  --counts2 "$COUNTS2" \
  --pos1 "$ALIGNED_POS" \
  --spatial2 "$SPATIAL2" \
  --outdir "$STCOMPARE_OUTDIR" \
  --sample_aligned "$SAMPLE_ALIGNED" \
  --sample_reference "$SAMPLE_REFERENCE"

deactivate_env

# final output summary
echo "Complete :D"
echo "Project directory:    $PROJECT_DIR"
echo "Tool directory:       $SCRIPT_DIR"
echo "Run directory:        $RUN_DIR"
echo "Landmarks:            $LANDMARK_DIR"
echo "Aligned coordinates:  $ALIGNED_POS"
echo "STalign outputs:      $STALIGN_OUTDIR"
echo "STcompare outputs:    $STCOMPARE_OUTDIR"