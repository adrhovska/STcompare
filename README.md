# STworkflow README
STworkflow.sh is a bash script which runs the full workflow of comparison of spatial transcriptomic data. It encompasses the LandmarkPicker, STalign and STcompare modules in sequential order. 

This workflow is designed for comparison of two 10x Visium ST tissue sections. It first aligns the source to the reference tissue with manually selected H&E landmarks (LandmarkPicker.py and STalignCode.py) and then continues with downstream spatial comparison with STcompare.R. 

The STalign parameters should be optimised according to the input samples individually to achieve the best fit. 

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

### STcompare
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
- No project specific paths are hardcoded in the script. Input data paths are passed through command-line arguments, and outputs are written to cwd or  selected project directory.

-------------------------------------------------------------------------------------------------------------------------------
The workflow file encompasses the following modules:

# LandmarkPicker
This Python script is used to manually select matching landmark pairs between two H&E stained images from 10x Visium tissue sections. They are saved as CSV files (y, x) and are designed to be further used for affine and LDDMM transformation in `STalignCode.py` in the pipeline. 

## Description
This code is a command line Python script that opens two H&E images, one source tissue image and one reference tissue image, and allows the user to click matching landmark points between them. For each landmark pair, the script first displays the source image and asks the user to click a landmark. It then closes it and displays the reference image asking the user to click the corresponding matching landmark. This is repeated for the prior selected number of landmark pairs (recommended 6 to 10).
The landmarks should be chosen on morphological basis, favourably equivalent anatomical features, or alternatively signatures such as indentations, tissue corners, folds, and internal cavities (lumen). They should ideally be spread out throughout the whole tissue section and not clustered in one area. 

The output files are saved in a directory named after the two samples and can be passed directly to `STalignCode.py` using the `--points1` and `--points2` arguments.

## Dependencies 
Python 3.9 and higher

The following Python packages: `pandas, matplotlib (pyplot and image), argparse, pathlib`

`LandmarkPicker.py` is a self-contained script and can be run by downloading it and passing arguments through the command line, replacing the example paths with file locations:

```text
python LandmarkPicker.py \
  --image1 /path/to/source_sample/spatial/tissue_hires_image.png \
  --image2 /path/to/reference_sample/spatial/tissue_hires_image.png \
  --sample_aligned SourceSample \
  --sample_reference ReferenceSample
```

After running the scripts asks the user to select the number of landmark pairs (6 to 10 is recommended), user should input the desired number of pairs. Afterwards the source image opens and the user is left to select the first of the pair, with the reference opening afterwards following first image closure, allowing the user to select the corresponding landmark. Repeat until all landmark pairs are selected. 

the code thus requires four arguments, those being: 
```bash
--image1 and 2 (.png source and reference hires images),
--sample_aligned (name of the source sample being aligned),
--sample_reference (name of the reference sample)
```

optional arguments include:
```bash
--project_dir (project directory where the landmark output folder will be created if working directory is wished to be overridden by the user)
```
## Outputs 
Outputs are saved in a directory named after the two samples in the project folder (has to be cd'). The directory contains source landmark CSV, reference landmark CSV and combined landmark CSV for QC. For `STalignCode.py` only the first two are essential, the last one serves as a check whether landmark order has been preserved (whether landmakr 1 in first tissue corresponds to landmark 1 in secod tissue).

## Additional Notes and QC
```text
-- matplotlib records clicked coordinates as x,y, but this script saves them as y,x because that is the coordinate format expected by STalignCode.py
-- in STalignCode.py would be passed as arguments pos1 (source tissue CSV) and pos2 (reference tissue CSV)
```
-------------------------------------------------------------------------------------------------------------------------------
The following text is expanded description and additional information about each of the modules in the pipeline. 

# STalign
This Python script aligns one 10x Visium tissue to a reference tissue based on manually selected H&E landmark pairs using `LandmarkPicker.py`. This fitted affine transform is then applied to Visium spot coordinates. 

## Description
This code is a command line Python script which aligns two ST samples from 10x Visium in order for them to be usable for STcompare. Using `LandmarkPicker.py` (separate Python script), landmarks are used to fit an affine transform from the source to the reference sample. This transform is then applied to the spot coordinates and generate an aligned coordinate file in the same coordinate space, which can then further be used for downstream spatial comparison using `STcompare`.

-- affine transformation: a geometric transformation that maps one coordinate system onto another while preserving straight lines. It includes translation, rotation, scaling, and shearing, but not local warping or bending. In this script, the affine transformation is fitted from manually selected landmark pairs from LandmarkPicker.py (has to be run separately) and then applied to all source Visium spot coordinates, aligning the source tissue section to the reference tissue section. 

## Dependencies
Python 3.9

The following Python packages: `numpy, pandas, matplotlib, argparse, json`

STalign is a self-contained script therefore can be run simply by downloading and passing on arguments using the command line in the following manner by replacing the paths with file locations: 

```text
python STalignCode.py \ 
    --pos1 /path/to/source_sample/spatial/tissue_positions.csv \ 
    --pos2 /path/to/reference_sample/spatial/tissue_positions.csv \ 
    --scale1 /path/to/source_sample/spatial/scalefactors_json.json \ 
    --scale2 /path/to/reference_sample/spatial/scalefactors_json.json \ 
    --points1 /path/to/landmarks/source_points.csv \ --> from LandmarkPicker
    --points2 /path/to/landmarks/reference_points.csv \ --> from LandmarkPicker
    --sample_aligned SourceSample \ 
    --sample_reference ReferenceSample 
```

The landmark files must contain matching landmark pairs in the same order. This means row 1 in points1 must correspond to row 1 in points2, row 2 to row 2 and so on. If the `LandmarkPicker.py` is used, such should be ensured. The expected format is CSV. 

```text
Example landmark CSV format:
y,x
320.5,510.2
480.1,1300.7
900.4,1500.3
1450.2,1200.8
1600.6,650.1
850.9,300.
```

the code thus requires **eight arguments**, those being: pos1 (CSV source tissue positions), pos2 (CSV reference tissue positions), scale1 (.json source scalefactors), scale2 (.json reference scalefactors), points1 and points2 (source and reference CSV y and x coordinates), sample_aligned (name of source tissue), sample_reference (name of reference tissue)

**optional arguments include**: project_dir (main project directory where output folders will be created), outdir (optional custom output directory), outname (optional custom name for the aligned barcode CSV)

## LDDMM parameter guide
In order for the non-linear fitting to run smoothly the parameters have to be optimised for the samples being compared. The STalign's LDDMM alignment is highly sensitive to the balance between several sigma parameters controlling the matching forces in the process, as well as to the step-size/regularisation parameters controlling how much the deformation can move. If these are disbalanced, not only will the alignment be incorrect but it might also fail entirely (empty output, crashed diagnostic QC plots, etc.)

'- sigmaM: noise standard deviation for the image matching term | Default: 1.0 | Smaller values mean a stronger pull toward matching image intensities between source and target.
'- sigmaB: noise standard deviation for the background term | Default: 2.0 | The same principle as SigmaM but for bcgrnd, made larger due to the matching of bcgrnd being not as important as that of the tissue.
'- sigmaA: noise standard deviation for the artifact term | Default: 5.0 | Artefacts have the weakest pull on the fit. Should be kept the largest.
'- sigmaP: noise standard deviation for landmark point-matching | Default: 20.0 | Controls how strongly the manually clicked landmarks are enforced. Making this too small (e.g. under 1) causes landmarks to pull far too aggressively relative to everything else, which can make the optimisation unstable. 
'- epV: step size for velocity field updates during optimization | Default: 1.0 | If you see runaway or divergent behavior (NaNs, blown-up deformation), this is usually the first thing to lower to around 1.
'- niter: total number of optimization iterations | Default: 1000
'- diffeo_start: the iteration at which LDDMM deformation begins, after an initial affine-only warm-up period | Default: 100.

**Note on parameter choices:** the sigma values act roughly as `1 / sigma²` in the
optimization's loss function so tightening `sigmaP` by 100x doesn't makes the landmark-matching **10,000x** stronger.
Therefore beware, that small changes to these values have an outsized,
nonlinear effect. Therefore adjustments in small jumps are advised.

If you run into an error, try troubleshooting:
1. Reset sigmas toward library defaults (`sigmaM=1, sigmaB=2, sigmaA=5, sigmaP=20`).
2. Lower `epV`.(around 1 is ideal)
3. Check landmark placement, such as landmarks that are too clustered, too few
   (fewer than ~6), or clicked in a mismatched order between source and
   reference images can poor initial affine which than fails LDDMM fitting

**Note on landmark reproducibility:** landmarks are chosen interactively
each run via `LandmarkPicker.py`, so results can vary slightly run to
run even with identical parameters. Therefore, before touching the parameters
it might be helpful to checl QC4 landmark QC first.

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

```text
Rscript STcompare.R \
  --counts1  /path/to/aligned_sample/filtered_feature_bc_matrix.h5 \
  --counts2  /path/to/reference_sample/filtered_feature_bc_matrix.h5 \
  --pos1     /path/to/aligned_sample(has to include barcodes).csv \
  --spatial2 /path/to/reference_sample/spatial \
  --outdir   /path/to/output_folder \
  --sample_aligned   (insert name) \
  --sample_reference (insert name)
```

the code thus requires **four arguments**, those being:
```bash
--counts1 (.h5 count matrix from Cell Ranger for the aligned sample),
--counts2 (.h5 count matrix from Cell Ranger for the reference sample),
--pos1 (CSV with columns barcode, x, y for the aligned sample),
--spatial2 (spatial/ folder from Cell Ranger for the reference sample)
```

**optional arguments include:**
```bash
--outdir (saving outputs)
--sample_aligned (name of aligned sample)
--sample_reference (name of reference sample)
--scale (hires or lowres img resolution)
--res (raster grid resolution, higher for more coarse and lower for finer)
--threads (CPU cores used)
```

## Outputs
Outputs are saved in a directory named after the two samples with subdirectories for results with correlation coefficients and p-values per gene tested, overlap of coordinates to check for scaling and correct alignment, raster graphs, correlation graphs (scatter of matched pixel expression values), regression graphs, and spatial map with color scheme based on predominant expression of a gene in sample above threshold (threshold is mean gene expression in the sample)

---

## Troubleshooting

### Missing file
Check that the source and reference folders contain the expected 10x Visium files:

```text
filtered_feature_bc_matrix.h5
spatial/tissue_hires_image.png
spatial/tissue_positions.csv
spatial/scalefactors_json.json
```

### Landmark files not found
Check whether `LandmarkPicker.py` saves landmark folders using the expected naming pattern:

```text
SampleAligned_paired_to_SampleReference
```

The workflow must point to the same landmark folder name.

### Conda activation error

Make sure that conda is installed and available in the terminal.

Also check that the environment names match those passed through or change accordingly in the argument passed:

```bash
--py_env
--r_env
```

### No matching barcodes

Check that the count matrix and position files come from the same sample.

Barcode suffix differences such as `-1` are handled in `STcompare.R`, but completely mismatched files will fail.

### Few genes found in both samples

Check whether the genes listed in `genes_of_interest` are present in both count matrices.


## Limitations
This workflow currently performs global affine alignment. Affine alignment can account for translation, rotation, scaling, and shearing, but it cannot correct local nonlinear deformation between tissue sections. Will be updated accordingly shortly.

The quality of the final alignment depends strongly on the quality and distribution of the manually selected landmarks. Landmarks should be clearly identifiable in both samples and spread across the full tissue section.

The STcompare results are limited to the genes defined in `genes_of_interest` inside `STcompare.R` as based on the paper https://doi.org/10.1038/s41587-026-03043-1. To analyse different genes, edit the gene list before running the workflow.

---

## Reproducibility Info

Each run saves the manually selected landmarks, the fitted affine transform, the aligned coordinate CSV, and all STcompare outputs inside one pair-specific output folder. This allows the full comparison to be traced back to the exact landmarks and affine transform used for alignment.

If the same sample names and output directory are reused, existing outputs may be overwritten. To preserve previous runs, use a new `--project_dir` or rename the output folder before rerunning the workflow.
