# load required libraries
suppressPackageStartupMessages({
  library(STcompare)
  library(SpatialExperiment)
  library(SEraster)
  library(hdf5r)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(jsonlite)
  library(argparser, quietly = TRUE)
})

## Parsing command line arguments with argparser
p <- arg_parser("STcompare")
# adding command line arguments
p <- add_argument(p, "--counts1", help = "Path to the aligned sample counts file")
p <- add_argument(p, "--counts2", help = "Path to the reference sample counts file")
p <- add_argument(p, "--pos1", help = "Path to the aligned positions file")
p <- add_argument(p, "--spatial2", help = "Path to the reference spatial directory")
p <- add_argument(p, "--outdir", help = "Output directory", default = "./STcompare_out")
p <- add_argument(p, "--scale", help = "Scale type: highres | lowres", default = "hires")
p <- add_argument(p, "--res", help = "Raster resolution", default = 150L, type = "integer")
p <- add_argument(p, "--threads", help = "Number of threads", default = 4L, type = "integer")
p <- add_argument(p, "--sample_aligned", help = "Name of the aligned sample", default = "Sample_1")
p <- add_argument(p, "--sample_reference", help = "Name of the reference sample", default = "Sample_2")
# parsing command line arguments
argv <- parse_args(p)
# validating arguments
required <- c("counts1", "counts2", "pos1", "spatial2")
missing <- required[sapply(required, function(k) is.na(argv[[k]]))]
if (length(missing) > 0) {
  print(paste("ERROR: Missing required argument(s):", paste("--", missing, collapse = ",")))
  quit(status = 1)
}

# seting names and print passed arguments
sample_aligned_name <- argv$sample_aligned
sample_reference_name <- argv$sample_reference

print(paste("counts1           :", argv$counts1))
print(paste("counts2           :", argv$counts2))
print(paste("pos1              :", argv$pos1))
print(paste("spatial2          :", argv$spatial2))
print(paste("outdir            :", argv$outdir))
print(paste("scale             :", argv$scale))
print(paste("res               :", argv$res))
print(paste("threads           :", argv$threads))
print(paste("sample_aligned    :", argv$sample_aligned))
print(paste("sample_reference  :", argv$sample_reference))

# creating output directories
# 1. create main output directory
dir.create(argv$outdir, showWarnings = FALSE, recursive = TRUE)
# 2. create comparison-specific subdirectory
comparison_name <- paste0(sample_aligned_name, "_vs_", sample_reference_name, "_", argv$scale, "_res", argv$res)
dir_comparison <- file.path(argv$outdir, comparison_name)
dir.create(dir_comparison, showWarnings = FALSE, recursive = TRUE)
# 3. create further subdirectories
output_names <- list()
output_dirs <- c("Results", "Coordinate_QC", "Raster_Plots", "Correlation_Plots", "Linear_Regression", "Pixel_Class")
output_dirs <- setNames(
  file.path(dir_comparison, output_names),
  output_names
)
for (d in output_names) {
  dir.create(output_dirs[[d]], showWarnings = FALSE, recursive = TRUE)
}

# defining genes of interest (from Supplementary material Extended Data Figure 8 b))
# and unlisting to then allow assignment of tissue type
genes_of_interest <- list(
  epithelial_genes = (c("KRT4", "KRT5", "IVL")),
  smooth_muscle_genes = (c("SMTN", "CALD1", "CSRP1", "TAGLN")),
  skeletal_muscle_genes = (c("TNNC1", "TNNC2", "ACTC1", "MYH8"))
) # Have to change downstream
genes_flat <- unlist(genes_of_interest, use.names = FALSE)

## Reader of aligned positions and Visium positions
# reads a CSV file containing aligned or Visium positions and returns a data frame with x and y coordinates
# checks for required columns and ensures that the coordinates are numeric and finite
# @path: path to the CSV file containing aligned positions
# @sample_name: name of the sample for error messages
# @type: type of positions ("aligned" or "visium")
# @scale_type: scale type for Visium positions ("hires" or "lowres")

read_positions <- function(path, sample_name, type = "visium", scale_type = "hires") {
  if (type == "aligned") {
    pos <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    required_cols <- c("barcode", "x", "y")
    if (!all(required_cols %in% colnames(pos))) {
      stop(paste(sample_name, "aligned file missing expected columns:", paste(required_cols, collapse = ", ")))
    }
    coord <- data.frame(x = as.numeric(pos$x), y = as.numeric(pos$y), row.names = pos$barcode)
  } else if (type == "visium") {
    pos_path <- file.path(path, "tissue_positions.csv")
    scale_path <- file.path(path, "scalefactors_json.json")
    pos <- read.csv(pos_path, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    required_cols <- c("pxl_row_in_fullres", "pxl_col_in_fullres")
    if (!all(required_cols %in% colnames(pos))) {
      stop(paste(sample_name, "tissue_positions.csv missing expected columns:", paste(required_cols, collapse = ", ")))
    }
    scales <- fromJSON(scale_path)
    scale_factor <- if (scale_type == "hires") scales$tissue_hires_scalef else scales$tissue_lowres_scalef
    cat(sample_name, scale_type, "scale factor:", scale_factor, "\n")
    coord <- data.frame(
      x = as.numeric(pos$pxl_col_in_fullres) * scale_factor,
      y = as.numeric(pos$pxl_row_in_fullres) * scale_factor,
      row.names = rownames(pos)
    )
  } else {
    stop("type must be 'aligned' or 'visium'")
  }
  if (any(!is.finite(coord$x)) || any(!is.finite(coord$y))) {
    stop(paste(sample_name, "coordinates contain non-finite values."))
  }
  return(coord)
}

# checking of coordinate system (might be removed, here due to issues with STalign coordinates and printing warnings --> will retain until fix)
check_coordinate_system <- function(coords1, coords2, sample_aligned_name, sample_reference_name) {
  ranges <- data.frame(
    sample = c(sample_aligned_name, sample_reference_name), min_x = c(min(coords1[, "x"]), min(coords2[, "x"])),
    max_x = c(max(coords1[, "x"]), max(coords2[, "x"])), min_y = c(min(coords1[, "y"]), min(coords2[, "y"])),
    max_y = c(max(coords1[, "y"]), max(coords2[, "y"]))
  )

  ranges$width <- ranges$max_x - ranges$min_x
  ranges$height <- ranges$max_y - ranges$min_y

  width_ratio <- max(ranges$width) / min(ranges$width)
  height_ratio <- max(ranges$height) / min(ranges$height)

  x_overlap <- max(coords1[, "x"]) >= min(coords2[, "x"]) &&
    max(coords2[, "x"]) >= min(coords1[, "x"])
  y_overlap <- max(coords1[, "y"]) >= min(coords2[, "y"]) &&
    max(coords2[, "y"]) >= min(coords1[, "y"])
  if (!x_overlap || !y_overlap) {
    print(paste("Coordinate ranges do not overlap, samples are not in the same aligned coordinate space."))
  }
  if (width_ratio > 3 || height_ratio > 3) {
    print(paste("Coordinate ranges differ in scale, one sample may be fullres while the other is hires/lowres."))
  }
}

## Matcher of counts to positions
# matches the barcodes in the counts matrix to the barcodes in the positions data frame
# returns a list containing the matched counts matrix and the corresponding coordinates
# -1 suffixes are removed from barcodes for matching if exact matches are not found
# @counts: counts matrix with barcodes as column names
# @pos: positions data frame with barcodes as row names
# @sample_name: name of the sample for error messages

match_counts_to_positions <- function(counts, pos, sample_name) {
  exact_common <- intersect(colnames(counts), rownames(pos))
  print(paste(sample_name, "exact barcode matches:", length(exact_common)))
  if (length(exact_common) > 0) {
    counts_m <- counts[, exact_common, drop = FALSE]
    pos_m <- pos[exact_common, , drop = FALSE]
  } else {
    count_key <- sub("-1$", "", colnames(counts))
    pos_key <- sub("-1$", "", rownames(pos))
    common_key <- intersect(count_key, pos_key)
    print(paste(sample_name, "cleaned barcode matches:", length(common_key)))
    if (length(common_key) == 0) stop(paste("No matching barcodes for", sample_name))
    count_idx <- match(common_key, count_key)
    pos_idx <- match(common_key, pos_key)
    counts_m <- counts[, count_idx, drop = FALSE]
    pos_m <- pos[pos_idx, , drop = FALSE]
    rownames(pos_m) <- colnames(counts_m)
  }
  coords <- cbind(x = as.numeric(pos_m$x), y = as.numeric(pos_m$y))
  rownames(coords) <- rownames(pos_m)
  if (ncol(counts_m) == 0 || nrow(coords) == 0) stop(paste(sample_name, " has zero matched spots."))
  if (!all(is.finite(coords))) stop(paste(sample_name, " has non-finite coordinates."))
  print(paste(sample_name, "final matched spots:", ncol(counts_m)))
  list(counts = counts_m, coords = coords)
}

# reading counts and positions, checking coordinate systems, and matching counts to positions using the defined functions
counts1 <- Read10X_h5(argv$counts1)
counts2 <- Read10X_h5(argv$counts2)
pos1 <- read_positions(argv$pos1, sample_aligned_name, type = "aligned")
pos2 <- read_positions(argv$spatial2, sample_reference_name, type = "visium", scale_type = argv$scale)
check_coordinate_system(pos1, pos2, sample_aligned_name, sample_reference_name)
samples <- list(
  list(counts = counts1, pos = pos1, name = sample_aligned_name),
  list(counts = counts2, pos = pos2, name = sample_reference_name)
)
matched <- lapply(samples, function(s) {
  match_counts_to_positions(s$counts, s$pos, s$name)
})
names(matched) <- c(sample_aligned_name, sample_reference_name)

counts1_matched <- matched[[sample_aligned_name]]$counts
coords1 <- matched[[sample_aligned_name]]$coords
counts2_matched <- matched[[sample_reference_name]]$counts
coords2 <- matched[[sample_reference_name]]$coords

# filtering genes of interest to those present in both count matrices
genes_flat <- genes_flat[
  genes_flat %in% rownames(counts1_matched) &
    genes_flat %in% rownames(counts2_matched)
]
if (length(genes_flat) == 0) {
  stop("None of the target genes are present in both count matrices.")
}
print(paste("Genes found in both samples:", length(genes_flat)))

# filtering the count matrices to only include the genes of interest
counts1_matched <- counts1_matched[genes_flat, , drop = FALSE]
counts2_matched <- counts2_matched[genes_flat, , drop = FALSE]

# creating a data frame for plotting the overlap of coordinates between the two samples
df_coords <- rbind(
  data.frame(coords1, sample = sample_aligned_name),
  data.frame(coords2, sample = sample_reference_name)
)
# coordinate overlap check plot
p_overlap <- ggplot(df_coords, aes(x = x, y = y, colour = sample)) +
  geom_point(size = 0.4, alpha = 0.5) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Coordinate overlap check", x = "x coordinate", y = "y coordinate", colour = "Sample")
ggsave(file.path(output_dirs[["Coordinate_QC"]], "Coordinate_QC.png"), p_overlap, width = 7, height = 6, dpi = 200)

## Object building
# creating SpatialExperiment objects for each sample using the matched counts and coordinates
spe_list <- setNames(
  mapply(function(counts, coords) {
    SpatialExperiment(assays = list(counts = counts), spatialCoords = coords)
  }, list(counts1_matched, counts2_matched), list(coords1, coords2), SIMPLIFY = FALSE),
  c(sample_aligned_name, sample_reference_name)
)

## Rasterization
rastList <- rasterizeGeneExpression(spe_list, assay_name = "counts", resolution = argv$res, square = FALSE)

## STcompare
sc <- spatialCorrelationGeneExp(rastList, nThreads = argv$threads)
ss <- spatialSimilarity(rastList)
# spatial correlation results for genes of interest
genes_in_sc <- genes_flat[genes_flat %in% rownames(sc)]
results <- sc[genes_in_sc, c("correlationCoef", "pValuePermuteX", "pValuePermuteY"), drop = FALSE]
# adding empirical p-value and cell type annotation to the results
results$empirical_pval <- pmax(results$pValuePermuteX, results$pValuePermuteY)
results$cell_type <- sapply(rownames(results), function(g) {
  names(which(sapply(genes_of_interest, function(grp) g %in% grp)))
})
print(results)

# saving the results to a CSV file in the results directory
write.csv(results, file.path(output_dirs[["Results"]], "Results_Table.csv"), row.names = TRUE)

# defining the raster assay to use for plotting
assays1 <- assayNames(rastList[[sample_aligned_name]])
assays2 <- assayNames(rastList[[sample_reference_name]])
common_assays <- intersect(assays1, assays2)
if (length(common_assays) == 0) {
  stop("No shared assay names between rasterised samples.")
}
rast_assay <- if ("counts" %in% common_assays) {
  "counts"
} else {
  warning(paste("Using first shared assay:", common_assays[1]))
  common_assays[1]
}
# defining shared x and y limits for plotting based on the coordinates of both samples
all_x <- Filter(is.finite, c(coords1[, "x"], coords2[, "x"]))
all_y <- Filter(is.finite, c(coords1[, "y"], coords2[, "y"]))
spatial_pad <- argv$res * 2
shared_xlim <- c(min(all_x) - spatial_pad, max(all_x) + spatial_pad)
shared_ylim <- c(min(all_y) - spatial_pad, max(all_y) + spatial_pad)
coord_label <- paste(argv$scale, "scale")
expr_label <- "rasterised raw counts"

## Builder of shared gene limits for raster plots
# rasterised values may differ between samples, finding the shared limits for each gene across both samples needed to ensure consistent color scaling in the plots
# @rastList: list of rasterised SpatialExperiment objects for each sample
# @gene: gene name for which to find the shared limits
# @assay_name: name of the assay to use for extracting raster values
# @name1: name of the first sample
# @name2: name of the second sample

get_shared_gene_lims <- function(rastList, gene, assay_name, name1, name2) {
  vals <- Filter(is.finite, as.numeric(c(assay(rastList[[name1]], assay_name)[gene, ], assay(rastList[[name2]], assay_name)[gene, ])))
  if (length(vals) == 0) stop("No finite raster values found for gene ", gene)
  limits <- range(vals)
  if (limits[1] == limits[2]) limits <- limits + c(-0.5, 0.5)
  return(limits)
}
# Creator of single raster plot for a given gene and sample
# @rast: rasterised SpatialExperiment object for the sample
# @name: name of the sample
# @gene: gene name for which to create the raster plot
# @gene_limits: shared limits for the gene across both samples for consistent color scaling
# @rast_assay: name of the assay to use for extracting raster values
# @shared_xlim: shared x-axis limits for the plot
# @shared_ylim: shared y-axis limits for the plot
# @coord_label: axis label describing the coordinate space (e.g. "hires scale")
# @expr_label: legend label describing the expression values (e.g. "rasterised raw counts")

make_single_raster <- function(rast, name, gene, gene_limits, rast_assay, shared_xlim, shared_ylim, coord_label, expr_label) {
  plotRaster(rast, assay_name = rast_assay, feature_name = gene, plotTitle = paste(name, "-", gene)) + scale_fill_viridis_c(limits = gene_limits, oob = scales::squish, name = paste0(gene, "\n", expr_label)) + coord_sf(xlim = shared_xlim, ylim = shared_ylim, expand = FALSE, clip = "off") + labs(
    x = paste("x coordinate (", coord_label, ")"),
    y = paste("y coordinate (", coord_label, ")")
  ) + theme(plot.margin = margin(10, 10, 10, 10))
}
## Joiner of raster plots for a single gene across two samples
# creates a side-by-side patchwork plot with a shared fill scale for comparability
# @gene: gene to plot
# @rastList: named list of rasterised SpatialExperiment objects
# @rast_assay: assay name to extract expression values from
# @name1: name of the aligned sample
# @name2: name of the reference sample
# @shared_xlim: shared x-axis limits for spatial comparability
# @shared_ylim: shared y-axis limits for spatial comparability
# @coord_label: axis label describing the coordinate space (e.g. "hires scale")
# @expr_label: legend label describing the expression values (e.g. "rasterised raw counts")

make_raster_pair <- function(gene, rastList, rast_assay, name1, name2, shared_xlim, shared_ylim, coord_label, expr_label) {
  gene_limits <- get_shared_gene_lims(rastList, gene, rast_assay, name1, name2)
  plots <- lapply(c(name1, name2), function(name) {
    make_single_raster(rastList[[name]], name, gene, gene_limits, rast_assay, shared_xlim, shared_ylim, coord_label, expr_label)
  })
  plots[[1]] + plots[[2]] + plot_layout(guides = "collect") & theme(legend.position = "right")
}

# loop over all genes of interest and save plots for each analysis type
# raster and correlation plots are only saved if the gene is present in the correlation results
save_plot <- function(plot, path, width, height) {
  ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
}
for (gene in genes_flat) {
  if (gene %in% genes_in_sc) {
    save_plot(
      make_raster_pair(
        gene, rastList, rast_assay, sample_aligned_name, sample_reference_name,
        shared_xlim, shared_ylim, coord_label, expr_label
      ),
      file.path(output_dirs[["Raster_Plots"]], paste0(gene, "_Raster.png")),
      width = 11, height = 5.5
    )
    save_plot(plotCorrelationGeneExp(rastList, sc, gene),
      file.path(output_dirs[["Correlation_Plots"]], paste0(gene, "_Correlation.png")),
      width = 10, height = 5
    )
  }
  save_plot(linearRegression(input = ss, gene = gene),
    file.path(output_dirs[["Linear_Regression"]], paste0(gene, "_LinearRegression.png")),
    width = 10, height = 5
  )
  save_plot(pixelClass(input = ss, gene = gene), file.path(output_dirs[["Pixel_Class"]], paste0(gene, "_PixelClass.png")), width = 10, height = 5)
}
