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