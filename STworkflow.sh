#!/bin/bash
set -euo pipefail

## User settings
# conda setup
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

# conda environments
PY_ENV="stalign_clean"
R_ENV="r_env"

# project folder
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(pwd)"

# parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source_dir)
      SOURCE_DIR="$2"
      shift 2
      ;;
    --reference_dir)
      REFERENCE_DIR="$2"
      shift 2
      ;;
    --sample_aligned)
      SAMPLE_ALIGNED="$2"
      shift 2
      ;;
    --sample_reference)
      SAMPLE_REFERENCE="$2"
      shift 2
      ;;
    --project_dir)
      PROJECT_DIR="$2"
      shift 2
      ;;
    --script_dir)
      SCRIPT_DIR="$2"
      shift 2
      ;;
    --py_env)
      PY_ENV="$2"
      shift 2
      ;;
    --r_env)
      R_ENV="$2"
      shift 2
      ;;
    --counts1)
      COUNTS1="$2"
      shift 2
      ;;
    --counts2)
      COUNTS2="$2"
      shift 2
      ;;
    --spatial2)
      SPATIAL2="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# check required arguments
if [[ -z "${SOURCE_DIR:-}" ]]; then
  echo "Missing required argument: --source_dir"
  usage
  exit 1
fi

if [[ -z "${REFERENCE_DIR:-}" ]]; then
  echo "Missing required argument: --reference_dir"
  usage
  exit 1
fi

if [[ -z "${SAMPLE_ALIGNED:-}" ]]; then
  echo "Missing required argument: --sample_aligned"
  usage
  exit 1
fi

if [[ -z "${SAMPLE_REFERENCE:-}" ]]; then
  echo "Missing required argument: --sample_reference"
  usage
  exit 1
fi

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

# helper checks
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

# conda activation helpers
activate_env() {
  set +u
  conda activate "$1"
  set -u
}
deactivate_env() {
  set +u
  conda deactivate
  set -u
}

# check tool scripts
need_file "$LANDMARK_PICKER"
need_file "$STALIGN_SCRIPT"
need_file "$STCOMPARE_SCRIPT"

# check input files
need_dir "$SOURCE_DIR"
need_dir "$REFERENCE_DIR"

need_file "$SOURCE_IMAGE"
need_file "$REFERENCE_IMAGE"

need_file "$SOURCE_POS"
need_file "$REFERENCE_POS"

need_file "$SOURCE_SCALE"
need_file "$REFERENCE_SCALE"

need_file "$COUNTS1"
need_file "$COUNTS2"
need_dir "$SPATIAL2"

# create main run directory
mkdir -p "$RUN_DIR"

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
  --pos1 "$SOURCE_POS" \
  --pos2 "$REFERENCE_POS" \
  --scale1 "$SOURCE_SCALE" \
  --scale2 "$REFERENCE_SCALE" \
  --points1 "$POINTS1" \
  --points2 "$POINTS2" \
  --sample_aligned "$SAMPLE_ALIGNED" \
  --sample_reference "$SAMPLE_REFERENCE" \
  --project_dir "$RUN_DIR" \
  --outdir "$STALIGN_OUTDIR"
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

echo "Complete :D"
echo "Project directory:    $PROJECT_DIR"
echo "Tool directory:       $SCRIPT_DIR"
echo "Run directory:        $RUN_DIR"
echo "Landmarks:            $LANDMARK_DIR"
echo "Aligned coordinates:  $ALIGNED_POS"
echo "STalign outputs:      $STALIGN_OUTDIR"
echo "STcompare outputs:    $STCOMPARE_OUTDIR"