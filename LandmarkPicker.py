from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

image1 = Path("/Users/adrhovska/Desktop/STdata/Native_1_ST/spatial/tissue_hires_image.png")
image2 = Path("/Users/adrhovska/Desktop/STdata/Native_2_ST/spatial/tissue_hires_image.png")
out1 = Path("/Users/adrhovska/Desktop/STdata/STcompare_code/points1.csv")
out2 = Path("/Users/adrhovska/Desktop/STdata/STcompare_code/points2.csv")
img1 = mpimg.imread(image1)
img2 = mpimg.imread(image2)
points1 = []
points2 = []
n_points = int(input("How many landmark pairs? Recommended 6-10: "))

for i in range(n_points):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img1)
    ax.set_title(f"Sample1/source: click landmark {i + 1}")
    clicked = plt.ginput(1, timeout=0)[0]
    plt.close(fig)

    x1, y1 = clicked
    points1.append({"y": y1, "x": x1})

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img2)
    ax.set_title(f"Sample2/reference: click matching landmark {i + 1}")
    clicked = plt.ginput(1, timeout=0)[0]
    plt.close(fig)

    x2, y2 = clicked
    points2.append({"y": y2, "x": x2})

pd.DataFrame(points1).to_csv(out1, index=False)
pd.DataFrame(points2).to_csv(out2, index=False)

print(f"Saved: {out1}")
print(f"Saved: {out2}")