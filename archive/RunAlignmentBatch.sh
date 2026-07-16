#!/opt/homebrew/bin/bash

SAMPLES_FILE="samples.txt"
PROJECT_DIR="$(pwd)/STcompare_batch"
SCRIPT_DIR="/Users/adrhovska/Desktop/STdata/STcompare_code"

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
"
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