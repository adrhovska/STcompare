#// ImageSplitter.py can be used to split the hires tissue image
# and the corresponding tissue positions into separate organoid folders for further analysis (separates organoids)
# pathways are hardcoded and have to be overwritten

import numpy as np
import pandas as pd
from pathlib import Path
import json
import shutil
from PIL import Image

# read paths to the directories and files 
block_spatial_dir = Path("/Users/adrhovska/Desktop/OCT_BLOCK_4/outs/spatial")
original_hires_image = block_spatial_dir / "tissue_hires_image.png"
scalefactors_file = block_spatial_dir / "scalefactors_json.json"
split_positions_dir = Path("/Users/adrhovska/Desktop/Split_Positions_B4")
outdir = Path("/Users/adrhovska/Desktop/BLOCK4_split_for_alignment")
margin = 150 #TODO: maybe smaller?

# read image
img = Image.open(original_hires_image).convert("RGB")
img_arr = np.array(img)
height, width = img_arr.shape[:2]

# read hires scale factor
with open(scalefactors_file, "r") as f:
    scalefactors = json.load(f)
hires_scale = scalefactors["tissue_hires_scalef"]
outdir.mkdir(parents=True, exist_ok=True)

# loop through organoid position files split by BarcodeMatcher.py
for pos_file in split_positions_dir.glob("*_tissue_positions.csv"):
    organoid_id = pos_file.stem.replace("_tissue_positions", "")
    positions = pd.read_csv(pos_file)

# convert full-resolution coordinates to hires-image coordinates
    x = positions["pxl_col_in_fullres"].astype(float) * hires_scale
    y = positions["pxl_row_in_fullres"].astype(float) * hires_scale

    # bounding box around this organoid
    x_min = max(int(x.min()) - margin, 0)
    x_max = min(int(x.max()) + margin, width)
    y_min = max(int(y.min()) - margin, 0)
    y_max = min(int(y.max()) + margin, height)

    # make white image of same size and then copy the organoid onto it from the original image
    masked = np.ones_like(img_arr) * 255
    masked[y_min:y_max, x_min:x_max, :] = img_arr[y_min:y_max, x_min:x_max, :]

    # output folder
    sample_spatial_dir = outdir / organoid_id / "spatial"
    sample_spatial_dir.mkdir(parents=True, exist_ok=True)

    # save masked same-size image
    Image.fromarray(masked).save(sample_spatial_dir / "tissue_hires_image.png")

    # save corresponding tissue positions
    positions.to_csv(sample_spatial_dir / "tissue_positions.csv", index=False)

    # copy scalefactors
    shutil.copy2(scalefactors_file, sample_spatial_dir / "scalefactors_json.json")

print("Saved:", organoid_id)