#!/bin/bash
set -euo pipefail

## User settings
# conda setup
source "$(conda info --base)/etc/profile.d/conda.sh"
# conda environments
py_env="stalign_clean"
r_env="r_env"

# project folder
script_dir="$(cd "$(dirname "${bash_source[0]}")" && pwd)"
project_dir="$(pwd)"
skip_landmarks=0

# parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source_dir)
      source_dir="$2"
      shift 2
      ;;
    --reference_dir)
      reference_dir="$2"
      shift 2
      ;;
    --sample_aligned)
      sample_aligned="$2"
      shift 2
      ;;
    --sample_reference)
      sample_reference="$2"
      shift 2
      ;;
    --project_dir)
      project_dir="$2"
      shift 2
      ;;
    --script_dir)
      script_dir="$2"
      shift 2
      ;;
    --py_env)
      py_env="$2"
      shift 2
      ;;
    --r_env)
      r_env="$2"
      shift 2
      ;;
    --counts1)
      counts1="$2"
      shift 2
      ;;
    --counts2)
      counts2="$2"
      shift 2
      ;;
    --spatial2)
      spatial2="$2"
      shift 2
      ;;
    --skip_landmarks)
      skip_landmarks=1
      shift 1
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
if [[ -z "${source_dir:-}" ]]; then
  echo "Missing required argument: --source_dir"
  usage
  exit 1
fi

if [[ -z "${reference_dir:-}" ]]; then
  echo "Missing required argument: --reference_dir"
  usage
  exit 1
fi

if [[ -z "${sample_aligned:-}" ]]; then
  echo "Missing required argument: --sample_aligned"
  usage
  exit 1
fi

if [[ -z "${sample_reference:-}" ]]; then
  echo "Missing required argument: --sample_reference"
  usage
  exit 1
fi

# derived paths
landmark_picker="${script_dir}/LandmarkPicker.py"
stalign_script="${script_dir}/STalignCode.py"
stcompare_script="${script_dir}/STcompare.R"
source_image="${source_dir}/spatial/tissue_hires_image.png"
reference_image="${reference_dir}/spatial/tissue_hires_image.png"
source_pos="${source_dir}/spatial/tissue_positions.csv"
reference_pos="${reference_dir}/spatial/tissue_positions.csv"
source_scale="${source_dir}/spatial/scalefactors_json.json"
reference_scale="${reference_dir}/spatial/scalefactors_json.json"
counts1="${counts1:-${source_dir}/filtered_feature_bc_matrix.h5}"
counts2="${counts2:-${reference_dir}/filtered_feature_bc_matrix.h5}"
spatial2="${spatial2:-${reference_dir}/spatial}"
pair_name="${sample_aligned}_paired_to_${sample_reference}"
landmark_dir="${project_dir}/landmarks/${pair_name}"
stalign_outdir="${project_dir}/STalign_outputs/${pair_name}"
stcompare_outdir="${project_dir}/STcompare_outputs/${pair_name}"
points1="${landmark_dir}/${sample_aligned}_points.csv"
points2="${landmark_dir}/${sample_reference}_points.csv"
aligned_pos="${stalign_outdir}/${sample_aligned}_aligned_to_${sample_reference}_barcodes.csv"

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
## Workflow
## Workflow

# 1: LandmarkPicker.py
  echo "Step 1: Running LandmarkPicker.py"
  echo "Click matching landmarks when the images open."
  conda activate "$py_env"
  python "$landmark_picker" \
    --image1 "$source_image" \
    --image2 "$reference_image" \
    --sample_aligned "$sample_aligned" \
    --sample_reference "$sample_reference" \
    --project_dir "$project_dir"
  conda deactivate

# 2: STalignCode.py
  echo "Step 2: Running STalignCode.py"
  conda activate "$py_env"
  python "$stalign_script" \
    --pos1 "$source_pos" \
    --pos2 "$reference_pos" \
    --scale1 "$source_scale" \
    --scale2 "$reference_scale" \
    --points1 "$points1" \
    --points2 "$points2" \
    --sample_aligned "$sample_aligned" \
    --sample_reference "$sample_reference" \
    --project_dir "$project_dir"
  conda deactivate


# step 3: STcompare.R
echo "Step 3: Running STcompare.R"
  conda activate "$r_env"
  Rscript "$stcomapre_script" \
    --counts1 "$counts1" \
    --counts2 "$counts2" \
    --pos1 "$aligned_pos" \
    --spatial2 "$spatial2" \
    --outdir "$stcompare_outdir" \
    --sample_aligned "$sample_aligned" \
    --sample_reference "$sample_reference"
  conda deactivate


echo "Complete :D"
echo "Project directory:    $PROJECT_DIR"
echo "Tool directory:       $SCRIPT_DIR"
echo "Landmarks:            $LANDMARK_DIR"
echo "Aligned coordinates:  $ALIGNED_POS"
echo "STalign outputs:      $STALIGN_OUTDIR"
echo "STcompare outputs:    $STCOMPARE_OUTDIR"