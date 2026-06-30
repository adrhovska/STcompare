# STworkflow README
STworkflow.sh is a bash script which runs the full workflow of comparison of spatial transcriptomic data. It encompasses the LandmarkPicker, STalign and STcompare modules in sequential order. 

This workflow is designed for comparison of two 10x Visium ST tissue sections. It first aligns the source to the reference tissue with manually selected H&E landmarks (LandmarkPicker.py and STalign.py) and then continues with downstream spatial comparison with STcompare.R. 

## Expected input
Each sample should be a standard 10x Visium output folder containing:

```text
Sample_ST/
|--- filtered_feature_bc_matrix.h5
|___ spatial/
    |--- tissue_hires_image.png
    |--- tissue_positions.csv
    |___ scalefactors_json.json
```

The workflow assumes the following default files for the two tissues:

SOURCE_DIR/filtered_feature_bc_matrix.h5
SOURCE_DIR/spatial/tissue_hires_image.png
SOURCE_DIR/spatial/tissue_positions.csv
SOURCE_DIR/spatial/scalefactors_json.json

REFERENCE_DIR/filtered_feature_bc_matrix.h5
REFERENCE_DIR/spatial/tissue_hires_image.png
REFERENCE_DIR/spatial/tissue_positions.csv
REFERENCE_DIR/spatial/scalefactors_json.json

## Requirements
The workflow uses two separate conda environments:
```bash
python_env
r_env
```

By default:
- `LandmarkPicker.py` and `STalignCode.py` run in python_env
- `STcompare.R` runs in r_env
The default environment names can be changed using:
```bash
--py_env PYTHON_CONDA_ENV
--r_env R_CONDA_ENV
```
## Usage
Make the workflow executable once:
```bash
chmod +x STworkflow.sh
```
Run the workflow:
```bash 
./STworkflow.sh \
  --source_dir /path/to/source_folder \
  --reference_dir /path/to/reference_folder \
  --sample_aligned SourceSample \
  --sample_reference ReferenceSample
```
the code thus requires four arguments, those being:
```bash 
--source_dir (path to the source 10x Visium folder, this is the sample being aligned),
--reference_dir (path to the reference 10x Visium folder, this is the template of the alignment),
--sample_aligned (name of the aligned sample),
--sample_reference (name of the reference sample)
```
optional arguments include:
```bash
--project_dir (output project directory, default is current working directory),
--script_dir (directory containing LandmarkPicker.py, STalignCode.py and STcompare.R, default is directory containing STworkflow.sh),
--py_env (Python conda environment name, default is python_env),
--r_env (R conda environment, default is r_env),
--counts1 (custom source count matrix path, default is filtered_feature_bc_matrix.h5),
--counts2 (custom reference count matrix path, default is REFERENCE_DIR/filtered_feature_bc_matrix.h5), 
--spatial2 (custom reference spatial folder path, default is REFERENCE_DIR/spatial)
```
## Outputs
The workflow creates one pair-specific output folder with Landmarks, STalign and STcompare folders. 

### Landmarks
Includes source and reference sample points and paired to files in CSV. The two essential files are (the paired file is used for QC, more below in LandmarkPicker section):
```text
Native_1_points.csv
Native_2_points.csv
```
These are passed into `STalignCode.py` as:
```text
--points1
--points2
```
### STalign

Includes CSV file, plots before and after alignment as well as landmark fit and affine transformation file:
```text
Native_1_aligned_to_Native_2_barcodes.csv
Native_1_vs_Native_2_before_alignment.png
Native_1_aligned_to_Native_2_after_alignment.png
Native_1_to_Native_2_manual_landmark_fit.png
Native_1_to_Native_2_manual_affine_transform.npz
```
The most important file is:
```text
Native_1_aligned_to_Native_2_barcodes.csv
```
This contains the aligned source spot coordinates and is passed into `STcompare.R`.

### STcmpare
The outputs of this folder contain:
```text
Results/
Coordinate_QC/
Raster_Plots/
Correlation_Plots/
Linear_Regression/
Pixel_Class/
```
The main result table is saved as:
```text
STcompare/Results/Results_Table.csv. This is a file key for comparison of the samples and is accompanied by multiple plots. 
```
## Notes
### Coordinate system 
- Manual landmarks are selected on `tissue_hires_image.png`. Therefore, `STalignCode.py` converts 10x Visium full-resolution spot coordinates into the hires image coordinate system. This ensures that landmark coordinates and spot coordinates are in the same coordinate system.
- The workflow always runs `LandmarkPicker.py`, so each run creates a new set of landmarks.
- If the same sample names and output folder are used again, previous landmark and alignment outputs may be overwritten.
- No project specific pathways are hardcoded in the script. Input data paths are passed through command-line arguments, and outputs are written to cwd or  selected project directory.

-------------------------------------------------------------------------------------------------------------------------------
The workflow file encompasses the following modules:

# LandmarkPicker
This Python script is used to manually select matching landmark pairs between two H&E stained images from 10x Visium tissue sections. They are saved as CSV files (y, x) and are designed to be further used for affine transformation in STalignCode.py. 

## Description
This code is a command line Python script that opens two H&E images, one source tissue image and one reference tissue image, and allows the user to click matching landmark points between them. For each landmark pair, the script first displays the source image and asks the user to click a landmark. It then displays the reference image and asks the user to click the corresponding matching landmark. This is repeated for the prior selected number of landmark pairs (recommended 6 to 10).

The output files are saved in a directory named after the two samples and can be passed directly to STalignCode.py using the `--points1` and `--points2` arguments.

landmark reccs: landmark pairs should be matching clear anatomical or structural points of the tissues visible in both H&E images. Such might include tissue corners, adventitia on the outside of the section, folds, holes and indentations, internal structures (lumen). They should ideally be spread out throughout the whole tissue section and not clustered in one area. 

## Dependencies 
Python 3.9 and higher
The following Python packages: pandas, matplotlib (pyplot and image), argparse, pathlib
LandmarkPicker.py is a self-contained script and can be run by downloading it and passing arguments through the command line, replacing the example paths with file locations:

```text
python LandmarkPicker.py \
  --image1 /path/to/source_sample/spatial/tissue_hires_image.png \
  --image2 /path/to/reference_sample/spatial/tissue_hires_image.png \
  --sample_aligned SourceSample \
  --sample_reference ReferenceSample
```

After running the scripts asks the user to select the number of landmark pairs (6 to 10 is recommended), user should input the desired number of pairs. Afterwards the source image opens and the user is left to select the first of the pair, with the reference opening afterwards allowing the user to select the corresponding landmark. Repeat until all landmarks are selected. 

the code thus requires four arguments, those being: image1 and 2 (.png source and reference hires images), sample_aligned (name of the source sample being aligned), sample_reference (name of the reference sample)

optional arguments include: project_dir (project directory where the landmark output folder will be created if working directory is wished to be overriden by the user)

## Outputs 
Outputs are saved in a directory named after the two samples in the project folder (has to be cd'). The directory contains source landmark CSV, reference landmark CSV and combined landmark CSV for QC. For STalignCode.py only the first two are essential, the last one serves as a check whether landmark order has been preserved.

## Additional Notes and QC
-- matplotlib records clicked coordinates as x,y, but this script saves them as y,x because that is the coordinate format expected by STalignCode.py
-- for H&E landmark picking from tissue_hires_image.png, the resulting landmark coordinates should be used with: --coord_scale hires in STalignCode.py
-- in STalignCode.py would be passed as arguments pos1 (source tissue CSV) and pos2 (reference tissue CSV)

-------------------------------------------------------------------------------------------------------------------------------

# STalign
This Python script aligns one 10x Visium tissue to a reference tissue based on manually selected H&E landmark pairs using LandmarkPicker.py. This fitted affine transform is then applied to Visium spot coordinates. 

## Description
This code is a command line Python script which aligns two ST samples from 10x Visium in order for them to be usable for STcompare. Using LandmarkPicker.py (separate Python script), landmarks are used to fit an affine transform from the source to the reference sample. This transform is then applied to the spot coordinates and generate an aligned coordinate file in the same coordinate space, which can then further be used for downstream spatial comparison using STcompare.

-- affine transformation: a geometric transformation that maps one coordinate system onto another while preserving straight lines. It includes translation, rotation, scaling, and shearing, but not local warping or bending. In this script, the affine transformation is fitted from manually selected landmark pairs from LandmarkPicker.py (has to be run separately) and then applied to all source Visium spot coordinates, aligning the source tissue section to the reference tissue section. 

## Dependencies
Python 3.9
The following Python packages: numpy, pandas, matplotlib, argparse, json 
STalign is a self-contained script therefore can be run simply by downloading and passing on arguments using the command line in the following manner by replacing the paths with file locations: 

python STalignCode.py \ 
    --pos1 /path/to/source_sample/spatial/tissue_positions.csv \ 
    --pos2 /path/to/reference_sample/spatial/tissue_positions.csv \ 
    --scale1 /path/to/source_sample/spatial/scalefactors_json.json \ 
    --scale2 /path/to/reference_sample/spatial/scalefactors_json.json \ 
    --points1 /path/to/landmarks/source_points.csv \ --> from LandmarkPicker
    --points2 /path/to/landmarks/reference_points.csv \ --> from LandmarkPicker
    --sample_aligned SourceSample \ 
    --sample_reference ReferenceSample \ 
    --coord_scale hires

The landmark files must contain matching landmark pairs in the same order. This means row 1 in points1 must correspond to row 1 in points2, row 2 to row 2 and so on. If the LandmarkPicker.py is used, such should be ensured. The expected format is CSV. 

Example landmark CSV format:
y,x
320.5,510.2
480.1,1300.7
900.4,1500.3
1450.2,1200.8
1600.6,650.1
850.9,300.
--> For normal H&E landmark picking from tissue_hires_image.png, use: --coord_scale hires (always)

the code thus requires eight arguments, those being: pos1 (CSV source tissue positions), pos2 (CSV reference tissue positions), scale1 (.json source scalefactors), scale2 (.json reference scalefactors), points1 and points2 (source and reference CSV y and x coordinates), sample_aligned (name of source tissue), sample_reference (name of reference tissue)

optional arguments include: project_dir (main project directory where output folders will be created), outdir (optional custom output directory), outname (optional custom name for the aligned barcode CSV), coord_scale (coordinate system used for spot and landmark coordinates, options are hires, lowres, fullres, or um), spot_diameter_um (Visium spot diameter in micrometres; default is 55.0)

## Outputs
Outputs are saved in a directory named after the two samples and contains CSV of aligned barcodes, before and after-alignment plot, manual landmark fit plot, affine transformation file in .npz format. 

## QC
The script prints the fitted affine matrix and landmark residuals, including mean, median, and maximum residual error. It also prints affine scale singular values and the affine determinant. For QC, check whether the landmark residuals are reasonably small, the affine determinant is positive and the affine scale values are not extremely small or extremely large. 

-------------------------------------------------------------------------------------------------------------------------------

# STcompare
This RScript compares spatial gene expression patterns of two Visium tissue sections.

## Description
This code is a command line RScript taking two 10x Visium spatial transcriptomics samples and comparing the spatial pattern and expression magnitude of genes specific for three tissue types: epithelium, smooth and skeletal muscle. This comparison is useful, for example, for assessing similarity between native and engineered tissue samples. 

## Dependencies
R 4.2 or higher
The following R packages: Seurat, ggplot2, patchwork, jsonlite, scales and SpatialExperiment, SEraster, STcompare and SummarizedExperiment from BiocManager 
STcompare is a self-contained script therefore can be run simply by downloading and passing on arguments using the command line in the following manner by replacing the paths with file locations: 

Rscript STcompare.R \
  --counts1  /path/to/aligned_sample/filtered_feature_bc_matrix.h5 \
  --counts2  /path/to/reference_sample/filtered_feature_bc_matrix.h5 \
  --pos1     /path/to/aligned_sample(has to include barcodes).csv \
  --spatial2 /path/to/reference_sample/spatial \
  --outdir   /path/to/output_folder \
  --sample_aligned   (insert name) \
  --sample_reference (insert name)

the code thus requires four arguments, those being: counts1 (.h5 count matrix from Cell Ranger for the aligned sample), counts2 (.h5 count matrix from Cell Ranger for the reference sample), pos1 (CSV with columns barcode, x, y for the aligned sample), spatial2 (spatial/ folder from Cell Ranger for the reference sample)

optional arguments include: outdir (saving outputs), sample_aligned (name of aligned sample), sample_reference (name of reference sample), scale (hires or lowres img resolution), res (raster grid resolution, higer for more coarse and lower for fineer), threads (CPU cores used)

## Outputs
Outputs are saved in a directory named after the two samples with subdirectories for results with correlation coefficients and p-values per gene tested, overlap of coordinates to check for scaling and correct alignment, raster graphs, correlation graphs (scatter of matched pixel expression values), regression graphs, and spatial map with color scheme based on predominant expression of a gene in sample above treshold (treshold is mean gene expression in the sample)