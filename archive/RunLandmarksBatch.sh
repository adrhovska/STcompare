#!/opt/homebrew/bin/bash

SAMPLES_FILE="samples.txt"
PROJECT_DIR="$(pwd)/STcompare_batch"
SCRIPT_DIR="/Users/adrhovska/Desktop/STdata/STcompare_code"
PY_ENV="python_env"

LANDMARK_PICKER="${SCRIPT_DIR}/LandmarkPicker.py"

mapfile -t NAMES < <(cut -f1 "$SAMPLES_FILE")
mapfile -t DIRS  < <(cut -f2 "$SAMPLES_FILE")
N=${#NAMES[@]}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$PY_ENV"

for ((i=0; i<N; i++)); do
  for ((j=i+1; j<N; j++)); do
    SAMPLE_ALIGNED="${NAMES[$i]}"
    SAMPLE_REFERENCE="${NAMES[$j]}"
    ALIGN_PAIR_NAME="${SAMPLE_ALIGNED}_aligned_to_${SAMPLE_REFERENCE}"
    LANDMARK_PAIR_NAME="${SAMPLE_ALIGNED}_paired_to_${SAMPLE_REFERENCE}"
    RUN_DIR="${PROJECT_DIR}/STcompare_outputs/${ALIGN_PAIR_NAME}"
    LANDMARK_DIR="${RUN_DIR}/landmarks/${LANDMARK_PAIR_NAME}"

    POINTS1="${LANDMARK_DIR}/${SAMPLE_ALIGNED}_points.csv"
    POINTS2="${LANDMARK_DIR}/${SAMPLE_REFERENCE}_points.csv"

    # skip if already picked
    if [[ -f "$POINTS1" && -f "$POINTS2" ]]; then
      echo "Skipping already-picked pair: $ALIGN_PAIR_NAME"
      continue
    fi
    mkdir -p "$RUN_DIR"
    python "$LANDMARK_PICKER" \
      --image1 "${DIRS[$i]}/spatial/tissue_hires_image.png" \
      --image2 "${DIRS[$j]}/spatial/tissue_hires_image.png" \
      --sample_aligned "$SAMPLE_ALIGNED" \
      --sample_reference "$SAMPLE_REFERENCE" \
      --project_dir "$RUN_DIR"
  done
done

conda deactivate