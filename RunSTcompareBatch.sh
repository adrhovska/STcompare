#!/opt/homebrew/bin/bash

# RunSTcompareBatch.sh
# Combines building of sample lists, batch generation and their alignment 
#
# Usage:
#   ./RunSTcompareBatch.sh                      # run all three phases
#   ./RunSTcompareBatch.sh --only list          # just build samples.txt
#   ./RunSTcompareBatch.sh --only landmarks     # just pick landmarks
#   ./RunSTcompareBatch.sh --only alignment     # just run alignment

set -euo pipefail

ORGANOIDS_DIR="/Users/adrhovska/Desktop/STdata/Organoids/Organoids_Split"
SCRIPT_DIR="/Users/adrhovska/Desktop/STdata/STcompare_code"
SAMPLES_FILE="samples.txt"
PROJECT_DIR="$(pwd)/STcompare_batch"
PY_ENV="python_env"
LANDMARK_PICKER="${SCRIPT_DIR}/LandmarkPicker.py"

ONLY="all"
if [[ "${1:-}" == "--only" ]]; then
  ONLY="$2"
fi

# phase 1: build sample list
build_samples_list() {
  > "$SAMPLES_FILE"
  for block_dir in "$ORGANOIDS_DIR"/BLOCK*_split_for_alignment; do
    for sample_dir in "$block_dir"/*/; do
      sample_name=$(basename "$sample_dir")
      if [[ -d "${sample_dir}spatial" ]]; then
        printf "%s\t%s\n" "$sample_name" "${sample_dir%/}" >> "$SAMPLES_FILE"
      fi
    done
  done
  echo "Wrote $(wc -l < "$SAMPLES_FILE" | tr -d ' ') samples to $SAMPLES_FILE"
}

# phase 2: pick landmarks for every pair 
run_landmarks_batch() {
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
}

# phase 3: run alignment for every pair
run_alignment_batch() {
  mapfile -t NAMES < <(cut -f1 "$SAMPLES_FILE")
  mapfile -t DIRS  < <(cut -f2 "$SAMPLES_FILE")
  N=${#NAMES[@]}

  PAIR_LOG="$PROJECT_DIR/alignment_completed.txt"
  mkdir -p "$PROJECT_DIR"
  touch "$PAIR_LOG"

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

      if grep -qx "$ALIGN_PAIR_NAME" "$PAIR_LOG"; then
        echo "Skipping already-aligned pair: $ALIGN_PAIR_NAME"
        continue
      fi
      if [[ ! -f "$POINTS1" || ! -f "$POINTS2" ]]; then
        echo "WARNING: no landmarks for $ALIGN_PAIR_NAME, skipping (run phase 1 first)"
        continue
      fi

      "$SCRIPT_DIR/STworkflow.sh" \
        --source_dir "${DIRS[$i]}" \
        --reference_dir "${DIRS[$j]}" \
        --sample_aligned "$SAMPLE_ALIGNED" \
        --sample_reference "$SAMPLE_REFERENCE" \
        --project_dir "$PROJECT_DIR" \
        --script_dir "$SCRIPT_DIR" \
        --skip_landmark_picking

      echo "$ALIGN_PAIR_NAME" >> "$PAIR_LOG"
    done
  done

  echo "All alignments and comparisons complete :D"
}

# run (which ones)
case "$ONLY" in
  list)       build_samples_list ;;
  landmarks)  run_landmarks_batch ;;
  alignment)  run_alignment_batch ;;
  all)
    build_samples_list
    run_landmarks_batch
    run_alignment_batch
    ;;
  *)
    echo "Unknown --only value: $ONLY (expected list|landmarks|alignment)" >&2
    exit 1
    ;;
esac
