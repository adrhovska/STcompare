# load required libraries
from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
default_dir = Path.cwd() # default directory is the directory where the script is run but can be overriden

# argument parsing 
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image1", required=True, help="Source tissue_hires_image.png")
    parser.add_argument("--image2", required=True, help="Reference tissue_hires_image.png")
    parser.add_argument("--sample_aligned", required=True, help="Source sample name")
    parser.add_argument("--sample_reference", required=True, help="Reference sample name")
    parser.add_argument("--project_dir", default=default_dir, type=Path, help="Output project directory. Defaults to the current working directory.")
    return parser.parse_args()

# clicking one landmark point from one image
def click_landmark(img, title):
    # creating and displaying window for the img 
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img)
    ax.set_title(title = "Pick landmark")
    # matplotlib returns coordinates as x, y
    clicked = plt.ginput(1, timeout=0)
    if not clicked:
        raise RuntimeError("No landmark was selected.")
    plt.close(fig)
    # storing clicked coordinates, first x then y 
    x, y = clicked
    # returning y first because this is needed for further analysis like this
    return y, x

def main():
    # reading terminal arguments and converting image/project paths from text strings to objects 
    args = parse_args()
    image1 = Path(args.image1)
    image2 = Path(args.image2)
    project_dir = Path(args.project_dir)

    # creating alignment and project drectory names 
    pair_name = f"{args.sample_aligned}_aligned_to_{args.sample_reference}"
    outdir = project_dir / "landmarks" / pair_name
    outdir.mkdir(parents=True, exist_ok=True)

    # output CSVs
    out1 = outdir / f"{args.sample_aligned}_points.csv"
    out2 = outdir / f"{args.sample_reference}_points.csv"
    out_combined = outdir / f"{pair_name}_landmark_pairs_combined.csv"

    # reading imgs into Python
    img1 = mpimg.imread(image1)
    img2 = mpimg.imread(image2)
    # point storage
    points1 = []
    points2 = []

    # defining the number of landmarks, usually around 6 to 10
    n_points = int(input("Choose amount of landmark points, usually 6 to 10: "))
    if n_points < 3:
        raise ValueError("Need at least 3 landmark pairs for affine alignment.")
    # looping over landmark pairs 
    for i in range(n_points):
        # clicking landmark on source images
        y1, x1 = click_landmark(img1, f"{args.sample_aligned}/source: click landmark {i + 1}")
        # saving points as dictionary with y first
        points1.append(
            {
                "landmark": i + 1,
                "y": y1,
                "x": x1,
            }
        )
        # clicking matching landmark on reference image
        y2, x2 = click_landmark(img2, f"{args.sample_reference}/reference: click matching landmark {i + 1}")
        # saving matching reference point as dictionary with y first
        points2.append(
            {
                "landmark": i + 1,
                "y": y2,
                "x": x2,
            }
        )

    # converting into pandas dfs
    df1 = pd.DataFrame(points1)
    df2 = pd.DataFrame(points2)
    # saving as CSVs for further analysis 
    df1[["y", "x"]].to_csv(out1, index=False)
    df2[["y", "x"]].to_csv(out2, index=False)
    # creates one combined table with both source and reference landmarks for QC
    combined = pd.DataFrame(
        {
            "landmark": df1["landmark"],
            f"{args.sample_aligned}_y": df1["y"],
            f"{args.sample_aligned}_x": df1["x"],
            f"{args.sample_reference}_y": df2["y"],
            f"{args.sample_reference}_x": df2["x"],
        }
    )
    # saving and printing directories of saving 
    combined.to_csv(out_combined, index=False)
    print(f"Source points:    {out1}")
    print(f"Reference points: {out2}")
    print(f"Combined file:    {out_combined}")

if __name__ == "__main__":
    main()