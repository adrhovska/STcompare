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
  library(gridExtra)
  library(argparser, quietly = TRUE)
})

# avoid macOS's quartz bitmap backend for ggsave()'s PNG output - quartz
# routes through the Aqua/Cocoa app framework and can pop open a visible
# R.app "Console" window on every fresh Rscript invocation (i.e. once per
# pair in a batch run). Cairo renders headless with no GUI involved.
if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

# parsing command line arguments with argparser
p <- arg_parser("STcompare")

p <- add_argument(p, "--counts1", help = "Path to the aligned sample counts file")
p <- add_argument(p, "--counts2", help = "Path to the reference sample counts file")
p <- add_argument(p, "--pos1", help = "Path to the aligned positions file")
p <- add_argument(p, "--spatial2", help = "Path to the reference spatial directory")
p <- add_argument(p, "--pos2", help = "Path to sample 2 positions (aligned CSV or visium spatial dir)")
p <- add_argument(p, "--type2", help = "Type of pos2: 'aligned' or 'visium'", default = "visium")
p <- add_argument(p, "--type1", help = "Type of pos1: 'aligned' or 'visium'", default = "aligned")
p <- add_argument(p, "--outdir", help = "Output directory", default = "./STcompare_out")
p <- add_argument(p, "--scale", help = "Scale type: highres | lowres", default = "hires")
p <- add_argument(p, "--res", help = "Raster resolution", default = 20L, type = "integer")
p <- add_argument(p, "--threads", help = "Number of threads", default = 4L, type = "integer")
p <- add_argument(p, "--sample_aligned", help = "Name of the aligned sample", default = "Sample_1")
p <- add_argument(p, "--sample_reference", help = "Name of the reference sample", default = "Sample_2")
argv <- parse_args(p)

# validating arguments
required <- c("counts1", "counts2", "pos1", "pos2")
missing <- required[sapply(required, function(k) is.na(argv[[k]]))]
if (length(missing) > 0) {
  print(paste("ERROR: Missing required argument(s):", paste("--", missing, collapse = ",")))
  quit(status = 1)
}

# setting names
sample_aligned_name <- argv$sample_aligned
sample_reference_name <- argv$sample_reference

# creating output directories
# 1. create main output directory
dir.create(argv$outdir, showWarnings = FALSE, recursive = TRUE)

# 2. create comparison-specific subdirectory
dir_comparison <- argv$outdir

# 3. create further subdirectories
# NEW: added Cluster_* directories to hold the tissue-type-cluster-level
# versions of the raster / correlation / regression / pixel-class plots,
# alongside (not replacing) the existing per-gene outputs
output_names <- c(
  "Results", "Raster_Plots", "Correlation_Plots", "Linear_Regression", "Pixel_Class",
  "Cluster_Raster_Plots", "Cluster_Correlation_Plots", "Cluster_Linear_Regression", "Cluster_Pixel_Class"
)
output_dirs <- setNames(file.path(dir_comparison, output_names), output_names)
for (d in output_names) {
  dir.create(output_dirs[[d]], showWarnings = FALSE, recursive = TRUE)
}

# defining genes of interest
# and unlisting to then allow assignment of tissue type
genes_of_interest <- list(
  progenitor_genes = c("HES1", "SOX2", "VIM"),
  maturation_genes = c("TH", "DCX", "MAP2"),
  patterning_genes = c("EN2", "NKX22", "WNT5A")
)
genes_flat <- unlist(genes_of_interest, use.names = FALSE)

# reading aligned positions and Visium positions
#   reads a CSV file containing aligned or Visium positions and returns a data frame with x and y coordinates
#   checks for required columns and ensures that the coordinates are numeric and finite
#   @path: path to the CSV file containing aligned positions
#   @sample_name: name of the sample for error messages
#   @type: type of positions ("aligned" or "visium")
#   @scale_type: scale type for Visium positions ("hires" or "lowres")

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

# matching of counts to positions
#   matches the barcodes in the counts matrix to the barcodes in the positions data frame
#   returns a list containing the matched counts matrix and the corresponding coordinates
#   @counts: counts matrix with barcodes as column names
#   @pos: positions data frame with barcodes as row names
#   @sample_name: name of the sample for error messages

match_counts_to_positions <- function(counts, pos, sample_name) {
  exact_common <- intersect(colnames(counts), rownames(pos))
  print(paste(sample_name, "exact barcode matches:", length(exact_common)))
  if (length(exact_common) > 0) {
    counts_m <- counts[, exact_common, drop = FALSE]
    pos_m <- pos[exact_common, , drop = FALSE]
  } else {
    count_key <- sub("-1$", "", colnames(counts)) # delete
    pos_key <- sub("-1$", "", rownames(pos))
    common_key <- intersect(count_key, pos_key)
    print(paste(sample_name, "cleaned barcode matches:", length(common_key)))
    if (length(common_key) == 0) stop(paste("No matching barcodes for", sample_name))
    count_idx <- match(common_key, count_key) # subset the cleaned
    pos_idx <- match(common_key, pos_key)
    counts_m <- counts[, count_idx, drop = FALSE]
    pos_m <- pos[pos_idx, , drop = FALSE]
    rownames(pos_m) <- colnames(counts_m)
  }

  # building the coordinates matrix for the matched spots
  coords <- cbind(x = as.numeric(pos_m$x), y = as.numeric(pos_m$y))
  rownames(coords) <- rownames(pos_m)

  # QC
  if (ncol(counts_m) == 0 || nrow(coords) == 0) stop(paste(sample_name, " has zero matched spots."))
  if (!all(is.finite(coords))) stop(paste(sample_name, " has non-finite coordinates."))
  print(paste(sample_name, "final matched spots:", ncol(counts_m)))
  list(counts = counts_m, coords = coords)
}

# reading counts and positions, checking coordinate systems, and matching counts to positions using the defined functions
counts1 <- Read10X_h5(argv$counts1)
counts2 <- Read10X_h5(argv$counts2)
pos1 <- read_positions(argv$pos1, sample_aligned_name, type = argv$type1, scale_type = argv$scale)
pos2 <- read_positions(argv$pos2, sample_reference_name, type = argv$type2, scale_type = argv$scale)
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

# NEW: restrict tissue-type cluster membership to genes that survived the
# both-samples filter above, and define the cluster names used throughout
# the cluster-level analysis below (progenitor_genes / maturation_genes /
# patterning_genes)
genes_of_interest <- lapply(genes_of_interest, function(g) g[g %in% genes_flat])
genes_of_interest <- genes_of_interest[sapply(genes_of_interest, length) > 0]
cluster_names <- names(genes_of_interest)
if (length(cluster_names) == 0) {
  warning("No tissue-type clusters have any surviving genes; cluster-level analysis will be skipped.")
}

# capturing total per-spot library size (ALL genes) before subsetting to the
# marker panel, so normalization is against real sequencing depth rather
# than just the sum of these few marker genes (which is often zero)
total_counts1 <- Matrix::colSums(counts1_matched)
total_counts2 <- Matrix::colSums(counts2_matched)

# filtering the count matrices to only include the genes of interest
counts1_matched <- counts1_matched[genes_flat, , drop = FALSE]
counts2_matched <- counts2_matched[genes_flat, , drop = FALSE]

# appending the total as an extra row so SEraster aggregates it into each
# pixel exactly like the marker genes
counts1_matched <- rbind(counts1_matched, `__TOTAL__` = total_counts1)
counts2_matched <- rbind(counts2_matched, `__TOTAL__` = total_counts2)

## Object building

# creating SpatialExperiment objects for each sample using the matched counts and coordinates
spe_list <- setNames(
  mapply(function(counts, coords) {
    SpatialExperiment(assays = list(counts = counts), spatialCoords = coords)
  }, list(counts1_matched, counts2_matched), list(coords1, coords2), SIMPLIFY = FALSE),
  c(sample_aligned_name, sample_reference_name)
)

# rasterization
rastList <- rasterizeGeneExpression(spe_list, assay_name = "counts", resolution = argv$res, square = FALSE)

for (name in names(rastList)) {
  mat <- assay(rastList[[name]])
  totals <- mat["__TOTAL__", ]
  keep <- setdiff(rownames(mat), "__TOTAL__")
  mat <- mat[keep, , drop = FALSE]
  totals[totals == 0] <- NA  # guard against any genuinely empty pixel
  mat <- t(t(mat) / totals) * 1e6  # CPM-normalised per-pixel expression

  # keep the original rasterised object's structure (pixel geometry, colData,
  # etc. that SEraster attaches and plotRaster() needs later) by subsetting
  # rather than rebuilding a bare SpatialExperiment from scratch
  gene_obj <- rastList[[name]][keep, ]
  assay(gene_obj) <- mat

  # NEW: build one pseudo-gene per tissue-type cluster by summing the
  # CPM-normalised expression of its member genes at each pixel. This turns
  # each cluster (progenitor / maturation / patterning) into a single
  # composite signal that can be pushed through exactly the same
  # correlation / similarity / raster pipeline used for individual genes,
  # so it is directly comparable to the paper's tissue-type-level results.
  if (length(cluster_names) > 0) {
    cluster_mat <- do.call(rbind, lapply(cluster_names, function(cl) {
      genes_here <- genes_of_interest[[cl]]
      colSums(mat[genes_here, , drop = FALSE])
    }))
    rownames(cluster_mat) <- cluster_names

    # clone the object's structure (same pixel geometry/colData) by
    # subsetting rows off gene_obj rather than constructing anything from
    # scratch, then swap in the cluster names/values; rbind merges it back
    # with the per-gene object since both share identical colData (pixels)
    cluster_obj <- gene_obj[seq_len(length(cluster_names)), ]
    rownames(cluster_obj) <- cluster_names
    assay(cluster_obj) <- cluster_mat

    rastList[[name]] <- rbind(gene_obj, cluster_obj)
  } else {
    rastList[[name]] <- gene_obj
  }
}

## STcompare
# sc / ss are computed over every row currently in rastList's assay, which
# now includes both the individual marker genes and the cluster pseudo-genes
sc <- spatialCorrelationGeneExp(rastList, nThreads = argv$threads)
ss <- spatialSimilarity(rastList) # installed STcompare version doesn't expose a fold-change threshold arg; uses its own default

# saves one row per pair with an overall similarity number
percent_similarity <- ss$similarityTable$percentSimilarity[
  match(genes_flat, ss$similarityTable$gene)
]

overall_similarity <- data.frame(
  sample_aligned = sample_aligned_name,
  sample_reference = sample_reference_name,
  mean_percent_similarity = mean(percent_similarity, na.rm = TRUE),
  n_genes_evaluated = sum(!is.na(percent_similarity))
)

write.csv(
  overall_similarity,
  file.path(output_dirs[["Results"]], "Overall_Similarity.csv"),
  row.names = FALSE
)
# spatial correlation results for genes of interest
genes_in_sc <- genes_flat[genes_flat %in% rownames(sc)]
results <- sc[genes_in_sc, c("correlationCoef", "pValuePermuteX", "pValuePermuteY"), drop = FALSE]

# adding empirical p-value and cell type annotation to the results
results$empirical_pval <- pmax(results$pValuePermuteX, results$pValuePermuteY)
results$cell_type <- sapply(rownames(results), function(g) {
  names(which(sapply(genes_of_interest, function(grp) g %in% grp)))
})
print(results)

# per-category summary (this is the mean of the individual-gene correlations
# within each cluster - a softer summary; see Cluster_Level_Results.csv /
# Cluster_Overall_Similarity.csv below for the true combined-signal metric)
category_summary <- aggregate(correlationCoef ~ cell_type, data = results, FUN = function(x) c(mean = mean(x, na.rm = TRUE), n = length(x)))
print(category_summary)
category_summary <- data.frame(cell_type = category_summary$cell_type, mean_correlation = round(category_summary$correlationCoef[, "mean"], 3), n_genes = category_summary$correlationCoef[, "n"])
write.csv(category_summary, file.path(output_dirs[["Results"]], "Category_Summary.csv"), row.names = FALSE)
# saving the results to a CSV file in the results directory
write.csv(results, file.path(output_dirs[["Results"]], "Results_Table.csv"), row.names = TRUE)

## NEW: tissue-type cluster-level analysis
# unlike Category_Summary.csv (average of independently-computed per-gene
# correlations), this runs the correlation/similarity computation directly
# on the summed cluster signal, and is the direct analogue of the per-gene
# Results_Table.csv / Overall_Similarity.csv at the tissue-type-cluster
# level - the numbers needed to build a per-cluster version of the paper's
# all-organoids similarity heatmap (see build_tissue_cluster_heatmaps.R)
if (length(cluster_names) > 0) {
  cluster_percent_similarity <- ss$similarityTable$percentSimilarity[
    match(cluster_names, ss$similarityTable$gene)
  ]
  cluster_overall_similarity <- data.frame(
    sample_aligned = sample_aligned_name,
    sample_reference = sample_reference_name,
    cluster = cluster_names,
    mean_percent_similarity = cluster_percent_similarity,
    n_genes_in_cluster = sapply(genes_of_interest, length)
  )
  write.csv(
    cluster_overall_similarity,
    file.path(output_dirs[["Results"]], "Cluster_Overall_Similarity.csv"),
    row.names = FALSE
  )

  clusters_in_sc <- cluster_names[cluster_names %in% rownames(sc)]
  results_cluster <- sc[clusters_in_sc, c("correlationCoef", "pValuePermuteX", "pValuePermuteY"), drop = FALSE]
  results_cluster$empirical_pval <- pmax(results_cluster$pValuePermuteX, results_cluster$pValuePermuteY)
  print(results_cluster)
  write.csv(results_cluster, file.path(output_dirs[["Results"]], "Cluster_Level_Results.csv"), row.names = TRUE)
} else {
  results_cluster <- NULL
  clusters_in_sc <- character(0)
}

# defining the raster assay to use for plotting
assays1 <- assayNames(rastList[[sample_aligned_name]])
assays2 <- assayNames(rastList[[sample_reference_name]])
common_assays <- intersect(assays1, assays2)
if (length(common_assays) == 0) {
  stop("No shared assay names between rasterised samples.") # should have counts from objectbuilder (could be removed?)
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
spatial_pad <- argv$res * 2 # rasterisation resolution
shared_xlim <- c(min(all_x) - spatial_pad, max(all_x) + spatial_pad)
shared_ylim <- c(min(all_y) - spatial_pad, max(all_y) + spatial_pad)
coord_label <- paste(argv$scale, "scale")
expr_label <- "rasterised raw counts"

# building shared gene limits for raster plots
#   rasterised values may differ between samples, finding the shared limits for each gene across both samples needed to ensure consistent color scaling in the plots
#   @rastList: list of rasterised SpatialExperiment objects for each sample
#   @gene: gene name for which to find the shared limits
#   @assay_name: name of the assay to use for extracting raster values
#   @name1: name of the first sample
#   @name2: name of the second sample

get_shared_gene_lims <- function(rastList, gene, assay_name, name1, name2) {
  vals <- Filter(is.finite, as.numeric(c(assay(rastList[[name1]], assay_name)[gene, ], assay(rastList[[name2]], assay_name)[gene, ])))
  if (length(vals) == 0) stop("No finite raster values found for gene ", gene)
  limits <- range(vals)
  if (limits[1] == limits[2]) limits <- limits + c(-0.5, 0.5)
  return(limits)
}

# creating single raster plot for a given gene and sample
#   @rast: rasterised SpatialExperiment object for the sample
#   @name: name of the sample
#   @gene: gene name for which to create the raster plot
#   @gene_limits: shared limits for the gene across both samples for consistent color scaling
#   @rast_assay: name of the assay to use for extracting raster values
#   @shared_xlim: shared x-axis limits for the plot
#   @shared_ylim: shared y-axis limits for the plot
#   @coord_label: axis label describing the coordinate space (e.g. "hires scale")
#   @expr_label: legend label describing the expression values (e.g. "rasterised raw counts")

make_single_raster <- function(rast, name, gene, gene_limits, rast_assay, shared_xlim, shared_ylim, coord_label, expr_label) {
  plotRaster(rast, assay_name = rast_assay, feature_name = gene, plotTitle = paste(name, "-", gene)) +
    scale_fill_viridis_c(limits = gene_limits, oob = scales::squish, name = paste0(gene, "\n", expr_label)) +
    coord_sf(xlim = shared_xlim, ylim = shared_ylim, expand = FALSE, clip = "off") + labs(
      x = paste("x coordinate (", coord_label, ")"),
      y = paste("y coordinate (", coord_label, ")")
    ) + theme(plot.margin = margin(10, 10, 10, 10))
}

# joining raster plots for a single gene across two samples
#   creates a side-by-side patchwork plot with a shared fill scale for comparability
#   @gene: gene to plot
#   @rastList: named list of rasterised SpatialExperiment objects
#   @rast_assay: assay name to extract expression values from
#   @name1: name of the aligned sample
#   @name2: name of the reference sample
#   @shared_xlim: shared x-axis limits for spatial comparability
#   @shared_ylim: shared y-axis limits for spatial comparability
#   @coord_label: axis label describing the coordinate space (e.g. "hires scale")
#   @expr_label: legend label describing the expression values (e.g. "rasterised raw counts")

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

## NEW: same set of diagnostic plots as above, but run on the tissue-type
# cluster pseudo-genes instead of individual marker genes
for (cl in cluster_names) {
  if (cl %in% clusters_in_sc) {
    save_plot(
      make_raster_pair(
        cl, rastList, rast_assay, sample_aligned_name, sample_reference_name,
        shared_xlim, shared_ylim, coord_label, expr_label
      ),
      file.path(output_dirs[["Cluster_Raster_Plots"]], paste0(cl, "_Raster.png")),
      width = 11, height = 5.5
    )
    save_plot(plotCorrelationGeneExp(rastList, sc, cl),
      file.path(output_dirs[["Cluster_Correlation_Plots"]], paste0(cl, "_Correlation.png")),
      width = 10, height = 5
    )
  }
  save_plot(linearRegression(input = ss, gene = cl),
    file.path(output_dirs[["Cluster_Linear_Regression"]], paste0(cl, "_LinearRegression.png")),
    width = 10, height = 5
  )
  save_plot(pixelClass(input = ss, gene = cl), file.path(output_dirs[["Cluster_Pixel_Class"]], paste0(cl, "_PixelClass.png")), width = 10, height = 5)
}

# final summary report
mean_corr <- mean(results$correlationCoef, na.rm = TRUE)
median_corr <- median(results$correlationCoef, na.rm = TRUE)
n_sig <- sum(results$empirical_pval < 0.05, na.rm = TRUE)
best_gene <- rownames(results)[which.max(results$correlationCoef)]
worst_gene <- rownames(results)[which.min(results$correlationCoef)]

fit_quality <- if (mean_corr >= 0.6) {
  "Strong spatial correspondence between samples for these markers."
} else if (mean_corr >= 0.3) {
  "Moderate spatial correspondence: some markers align well, others may need review."
} else {
  "Weak spatial correspondence (try reviewing alignment quality first, check QC outputs from STalign)"
}

# NEW: tissue-type cluster line for the summary text
cluster_summary_text <- if (!is.null(results_cluster) && nrow(results_cluster) > 0) {
  paste0(
    "Tissue-type cluster correlations:\n",
    paste(sprintf("  %s: r=%.3f (p=%.3f)", rownames(results_cluster), results_cluster$correlationCoef, results_cluster$empirical_pval), collapse = "\n")
  )
} else {
  "Tissue-type cluster correlations: none computed"
}

# key results in gene expression matching
summary_text <- paste(
  paste0(sample_aligned_name, " vs ", sample_reference_name),
  paste0("Matched spots: ", ncol(counts1_matched), " / ", ncol(counts2_matched)),
  paste0("Genes analyzed: ", nrow(results), "   Significant (p<0.05): ", n_sig, " / ", nrow(results)),
  paste0("Mean r = ", round(mean_corr, 3), "   Median r = ", round(median_corr, 3)),
  paste0(
    "Best: ", best_gene, " (r=", round(max(results$correlationCoef, na.rm = TRUE), 3),
    ")   Worst: ", worst_gene, " (r=", round(min(results$correlationCoef, na.rm = TRUE), 3), ")"
  ),
  "",
  cluster_summary_text,
  "",
  fit_quality,
  sep = "\n"
)

text_panel <- ggplot() +
  annotate("text", x = 0, y = 0, label = summary_text, hjust = 0, vjust = 1, size = 4, family = "mono") +
  xlim(0, 10) +
  ylim(-8, 1) +
  theme_void()

# results table
table_panel <- tableGrob(
  round(results[, c("correlationCoef", "empirical_pval")], 3),
  theme = ttheme_minimal(base_size = 8)
)

# reusing the best-correlated gene's raster pair and correlation plot
best_gene_limits <- get_shared_gene_lims(rastList, best_gene, rast_assay, sample_aligned_name, sample_reference_name)
key_raster <- make_raster_pair(
  best_gene, rastList, rast_assay, sample_aligned_name, sample_reference_name,
  shared_xlim, shared_ylim, coord_label, expr_label
)
key_corr_plot <- plotCorrelationGeneExp(rastList, sc, best_gene)

# combination into final document and directing it to the right saving point
final_panel <- (text_panel | table_panel) / key_raster / key_corr_plot +
  plot_layout(heights = c(1, 1.3, 1)) +
  plot_annotation(title = paste("STcompare Summary:", sample_aligned_name, "aligned to", sample_reference_name))

summary_path <- file.path(output_dirs[["Results"]], "Summary_Report.pdf")
ggsave(summary_path, final_panel, width = 12, height = 14, dpi = 300, bg = "white")

print(paste("Summary report saved to:", summary_path))