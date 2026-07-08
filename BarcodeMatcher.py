import pandas as pd
from pathlib import Path

# hardcode paths to the compared barcodes: BLOCK1
full_block_barcodes = Path("/Users/adrhovska/Desktop/OCT_BLOCK_4/outs/filtered_feature_bc_matrix/barcodes.tsv")
block_tissue_positions = Path("/Users/adrhovska/Desktop/OCT_BLOCK_4/outs/spatial/tissue_positions.csv")

# folder with the complete set of sorted organoid barcode files
sorted_barcode_folder = Path("/Users/adrhovska/Desktop/OCT_barcodes")

# output folder for split tissue-position files
outdir = Path("/Users/adrhovska/Desktop/Split_Positions")


def read_barcodes(path):
    path = Path(path)
    lines = path.read_text().splitlines()

    barcodes = set()

    for line in lines:
        if not line:
            continue
        barcodes.add(line)
    return barcodes

# create output folder
outdir.mkdir(parents=True, exist_ok=True)

# read the barcodes from the full block
whole_block_barcodes = read_barcodes(full_block_barcodes)

# read tissue positions from the full block
positions = pd.read_csv(block_tissue_positions)
positions["barcode"] = positions["barcode"].astype(str)

summary = []

# loop all of the sorted barcode files
for barcode_file in sorted_barcode_folder.glob("*barcodes*"):
    organoid_id = barcode_file.stem.replace("_barcodes", "")

    sorted_barcodes = read_barcodes(barcode_file)

    # barcodes shared between this organoid list and the whole block
    overlap = sorted_barcodes & whole_block_barcodes

    # require the whole barcode list to be present in this block
    if len(overlap) != len(sorted_barcodes):
        summary.append({
            "organoid_id": organoid_id,
            "barcode_file": barcode_file.name,
            "n_sorted_barcodes": len(sorted_barcodes),
            "n_overlap_with_block": len(overlap),
            "n_missing_from_block": len(sorted_barcodes) - len(overlap),
            "n_positions_saved": 0,
            "saved": False,
        })
        continue

    # split tissue_positions.csv
    subset = positions[positions["barcode"].isin(overlap)].copy()
    subset["organoid_id"] = organoid_id
    subset["barcode_source_file"] = barcode_file.name

    # save one CSV per organoid
    output_file = outdir / f"{organoid_id}_tissue_positions.csv"
    subset.to_csv(output_file, index=False)

    summary.append({
        "organoid_id": organoid_id,
        "barcode_file": barcode_file.name,
        "n_sorted_barcodes": len(sorted_barcodes),
    })

# save summary
summary = pd.DataFrame(summary)
summary.to_csv(outdir / "split_summary.csv", index=False)
