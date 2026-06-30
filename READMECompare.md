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