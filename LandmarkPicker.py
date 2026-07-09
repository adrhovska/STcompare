# load required libraries
from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
default_dir = Path.cwd() # default directory is the directory where the script is run but can be overriden 

#// argument parsing (use argparse library available online: https://docs.python.org/3/library/argparse.html)
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image1", required=True, help="Source tissue_hires_image.png")
    parser.add_argument("--image2", required=True, help="Reference tissue_hires_image.png")
    parser.add_argument("--sample_aligned", required=True, help="Source sample name")
    parser.add_argument("--sample_reference", required=True, help="Reference sample name")
    parser.add_argument("--project_dir", default=default_dir, type=Path, help="Output project directory. Defaults to the current working directory.") # type = Path helps it from returning a string to the path being stored as an object
    return parser.parse_args()

#// function to click landmarks in images
# creates a displaying window for the image and opens it
# waits for the user to click on the landmark and returns the coordinates of the clicked point
# returns coordinates in y, x order (STalign requires it in this order, y are rows and x columns)
# closes image and displays the next one 
# @param img: image to display
# @param title: title of the displayed image 
def click_landmark(img, title):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img)
    ax.set_title(title)
    clicked = plt.ginput(1, timeout=0)
    if not clicked:
        raise RuntimeError("No landmark was selected.")
    plt.close(fig)
    x, y = clicked[0]
    return y, x

# main body
def main():
    # parsing arguments and creating paths for images and project directory
    args = parse_args()
    image1 = Path(args.image1)
    image2 = Path(args.image2)
    project_dir = Path(args.project_dir)

    # creating landmarks output directory for the two samples
    outdir = project_dir / "landmarks" / f"{args.sample_aligned}_paired_to_{args.sample_reference}"
    outdir.mkdir(parents=True, exist_ok=True)

    # creating paths for the output CSVs files
    out1 = outdir / f"{args.sample_aligned}_points.csv"
    out2 = outdir / f"{args.sample_reference}_points.csv"
    out_combined = outdir / f"{args.sample_aligned}_and_{args.sample_reference}_landmark_pairs_combined.csv"

    # reading imgs into Python and storing them as arrays (for now empty lists)
    img1 = mpimg.imread(image1)
    img2 = mpimg.imread(image2)
    points1 = []
    points2 = []

    # setting the number of landmarks to choose on the image (3 is the minimum allowed due to affine alignment)
    n_points = int(input("Choose amount of landmark points, usually 6 to 10: "))
    if n_points < 3:
    # looping over the chosen landmarks and storing the coordinate outputs in the lists (saved as dictionaries) 
    # repeating the same for the second image 
        y1, x1 = click_landmark(img1, f"{args.sample_aligned}/source: click landmark {i + 1}")
        raise ValueError("Need at least 3 landmark pairs for affine alignment.")
    for i in range(n_points): 
        points1.append(
            {
                "landmark": i + 1,
                "y": y1,
                "x": x1,
            }
        )
        y2, x2 = click_landmark(img2, f"{args.sample_reference}/reference: click matching landmark {i + 1}")
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
    # creating one combined table with both source and reference landmarks for QC (check correct pairing of the landmarks)
    combined = pd.DataFrame(
        {
            "landmark": df1["landmark"], # can take from any df essentially 
            f"{args.sample_aligned}_y": df1["y"],
            f"{args.sample_aligned}_x": df1["x"],
            f"{args.sample_reference}_y": df2["y"],
            f"{args.sample_reference}_x": df2["x"],
        }
    )
    combined.to_csv(out_combined, index=False)

if __name__ == "__main__":
    main()