#!/opt/homebrew/bin/bash

ORGANOIDS_DIR="/Users/adrhovska/Desktop/STdata/Organoids/Organoids_Split"
OUT_FILE="samples.txt"

> "$OUT_FILE"

for block_dir in "$ORGANOIDS_DIR"/BLOCK*_split_for_alignment; do
  for sample_dir in "$block_dir"/*/; do
    sample_name=$(basename "$sample_dir")
    if [[ -d "${sample_dir}spatial" ]]; then
      printf "%s\t%s\n" "$sample_name" "${sample_dir%/}" >> "$OUT_FILE"
    fi
  done
done