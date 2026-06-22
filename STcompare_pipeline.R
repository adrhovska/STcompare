suppressPackageStartupMessages({
  library(STcompare)
  library(SpatialExperiment)
  library(SEraster)
  library(hdf5r)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    counts1 = NULL,
    counts2 = NULL,
    pos1    = NULL,
    spatial2 = NULL,
    outdir  = "./STcompare_out",
    scale   = "hires",
    res     = 150L,
    threads = 4L
  )

  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    val <- if (i + 1 <= length(args)) args[i + 1] else stop(paste("Missing value for --", key))
    defaults[[key]] <- val
    i <- i + 2
  }

  defaults$res     <- as.integer(defaults$res)
  defaults$threads <- as.integer(defaults$threads)

  required <- c("counts1", "counts2", "pos1", "spatial2")
  missing  <- required[sapply(required, function(k) is.null(defaults[[k]]))]
  if (length(missing)) {
    cat("ERROR: Missing required argument(s):", paste0("--", missing, collapse = ", "), "\n\n")
    cat("USAGE:\n")
    cat("  Rscript STcompare_pipeline.R \\\n")
    cat("    --counts1  <Native_1.h5> \\\n")
    cat("    --counts2  <Native_2.h5> \\\n")
    cat("    --pos1     <Native_1_aligned_WITH_BARCODES.csv> \\\n")
    cat("    --spatial2 <Native_2_ST/spatial/> \\\n")
    cat("    [--outdir  ./STcompare_out] \\\n")
    cat("    [--scale   hires|lowres] \\\n")
    cat("    [--res     150] \\\n")
    cat("    [--threads 4]\n")
    quit(status = 1)
  }

  defaults
}

cfg <- parse_args(args)

cat("\n=== STcompare Pipeline ===\n")
cat("counts1  :", cfg$counts1,  "\n")
cat("counts2  :", cfg$counts2,  "\n")
cat("pos1     :", cfg$pos1,     "\n")
cat("spatial2 :", cfg$spatial2, "\n")
cat("outdir   :", cfg$outdir,   "\n")
cat("scale    :", cfg$scale,    "\n")
cat("res      :", cfg$res,      "\n")
cat("threads  :", cfg$threads,  "\n\n")

dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)

epithelial_genes      <- c("KRT4", "KRT5", "IVL")
smooth_muscle_genes   <- c("SMTN", "CALD1", "CSRP1", "TAGLN")
skeletal_muscle_genes <- c("TNNC1", "TNNC2", "ACTC1", "MYH8")

genes_of_interest <- c(epithelial_genes, smooth_muscle_genes, skeletal_muscle_genes)

get_gene_expression <- function(counts) {
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) return(counts[["Gene Expression"]])
    return(counts[[1]])
  }
  counts
}

read_aligned_with_barcodes <- function(path) {
  pos <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("barcode", "x", "y") %in% colnames(pos)))
    stop("Aligned file must contain columns: barcode, x, y")
  pos_scaled <- data.frame(x = as.numeric(pos$x), y = as.numeric(pos$y), row.names = pos$barcode)
  if (any(!is.finite(pos_scaled$x)) || any(!is.finite(pos_scaled$y)))
    stop("Aligned Native 1 coordinates contain non-finite values.")
  pos_scaled
}

read_visium_positions_scaled <- function(spatial_dir, scale_type = "hires") {
  pos_path   <- file.path(spatial_dir, "tissue_positions.csv")
  scale_path <- file.path(spatial_dir, "scalefactors_json.json")
  pos        <- read.csv(pos_path, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(pos)))
    stop("tissue_positions.csv missing expected columns.")
  scales       <- jsonlite::fromJSON(scale_path)
  scale_factor <- if (scale_type == "hires") scales$tissue_hires_scalef else scales$tissue_lowres_scalef
  cat("Native 2", scale_type, "scale factor:", scale_factor, "\n")
  pos_scaled <- data.frame(
    x = as.numeric(pos$pxl_col_in_fullres) * scale_factor,
    y = as.numeric(pos$pxl_row_in_fullres) * scale_factor,
    row.names = rownames(pos)
  )
  if (any(!is.finite(pos_scaled$x)) || any(!is.finite(pos_scaled$y)))
    stop("Native 2 scaled coordinates contain non-finite values.")
  pos_scaled
}

match_counts_to_positions <- function(counts, pos, sample_name) {
  exact_common <- intersect(colnames(counts), rownames(pos))
  cat(sample_name, "exact barcode matches:", length(exact_common), "\n")
  if (length(exact_common) > 0) {
    counts_m <- counts[, exact_common, drop = FALSE]
    pos_m    <- pos[exact_common, , drop = FALSE]
  } else {
    count_key <- sub("-1$", "", colnames(counts))
    pos_key   <- sub("-1$", "", rownames(pos))
    common_key <- intersect(count_key, pos_key)
    cat(sample_name, "cleaned barcode matches:", length(common_key), "\n")
    if (length(common_key) == 0) stop(paste0("No matching barcodes for ", sample_name))
    count_idx <- match(common_key, count_key)
    pos_idx   <- match(common_key, pos_key)
    counts_m  <- counts[, count_idx, drop = FALSE]
    pos_m     <- pos[pos_idx, , drop = FALSE]
    rownames(pos_m) <- colnames(counts_m)
  }
  coords <- cbind(x = as.numeric(pos_m$x), y = as.numeric(pos_m$y))
  rownames(coords) <- rownames(pos_m)
  if (ncol(counts_m) == 0 || nrow(coords) == 0) stop(paste0(sample_name, " has zero matched spots."))
  if (!all(is.finite(coords))) stop(paste0(sample_name, " has non-finite coordinates."))
  cat(sample_name, "final matched spots:", ncol(counts_m), "\n")
  list(counts = counts_m, coords = coords)
}

cat("Loading count matrices...\n")
counts1 <- get_gene_expression(Read10X_h5(cfg$counts1))
counts2 <- get_gene_expression(Read10X_h5(cfg$counts2))
cat("Native 1 dims:", dim(counts1), "\n")
cat("Native 2 dims:", dim(counts2), "\n")

cat("\nLoading positions...\n")
pos1 <- read_aligned_with_barcodes(cfg$pos1)
pos2 <- read_visium_positions_scaled(cfg$spatial2, scale_type = cfg$scale)

cat("\nMatching barcodes to positions...\n")
matched1 <- match_counts_to_positions(counts1, pos1, "Native_1")
matched2 <- match_counts_to_positions(counts2, pos2, "Native_2")

counts1_matched <- matched1$counts
coords1         <- matched1$coords
counts2_matched <- matched2$counts
coords2         <- matched2$coords


genes_of_interest <- genes_of_interest[
  genes_of_interest %in% rownames(counts1_matched) &
  genes_of_interest %in% rownames(counts2_matched)
]

if (length(genes_of_interest) == 0)
  stop("None of the target genes are present in both count matrices.")

cat("\nGenes found in both samples:\n")
print(genes_of_interest)

counts1_matched <- counts1_matched[genes_of_interest, , drop = FALSE]
counts2_matched <- counts2_matched[genes_of_interest, , drop = FALSE]


df_coords <- rbind(
  data.frame(x = coords1[, "x"], y = coords1[, "y"], sample = "Native_1_aligned"),
  data.frame(x = coords2[, "x"], y = coords2[, "y"], sample = "Native_2_scaled")
)

p_overlap <- ggplot(df_coords, aes(x = x, y = y, colour = sample)) +
  geom_point(size = 0.4, alpha = 0.5) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Coordinate overlap check")

ggsave(file.path(cfg$outdir, "coordinate_overlap.png"), p_overlap, width = 7, height = 6, dpi = 200)
cat("Saved: coordinate_overlap.png\n")


cat("\nBuilding SpatialExperiment objects...\n")
spe1 <- SpatialExperiment(assays = list(counts = counts1_matched), spatialCoords = coords1)
spe2 <- SpatialExperiment(assays = list(counts = counts2_matched), spatialCoords = coords2)

cat("Rasterising gene expression (resolution =", cfg$res, ")...\n")
rastList <- SEraster::rasterizeGeneExpression(
  list(Native_1 = spe1, Native_2 = spe2),
  assay_name = "counts",
  resolution = cfg$res,
  square     = FALSE
)


cat("\nComputing spatial correlation...\n")
sc <- spatialCorrelationGeneExp(
  list(Native_1 = rastList$Native_1, Native_2 = rastList$Native_2),
  nThreads = cfg$threads
)

cat("Computing spatial similarity...\n")
ss <- spatialSimilarity(
  list(Native_1 = rastList$Native_1, Native_2 = rastList$Native_2)
)


genes_in_sc       <- genes_of_interest[genes_of_interest %in% rownames(sc)]
results           <- sc[genes_in_sc, c("correlationCoef", "pValuePermuteX", "pValuePermuteY"), drop = FALSE]
results$empirical_pval <- pmax(results$pValuePermuteX, results$pValuePermuteY)
results$cell_type <- ifelse(
  rownames(results) %in% epithelial_genes, "Epithelial",
  ifelse(rownames(results) %in% smooth_muscle_genes, "Smooth Muscle", "Skeletal Muscle")
)

cat("\n=== Results ===\n")
print(results)

write.csv(results, file.path(cfg$outdir, "results_table.csv"))
cat("Saved: results_table.csv\n")


cat("\nSaving plots...\n")

for (gene in genes_of_interest) {

  # Raster plots
  p1 <- plotRaster(rastList$Native_1, feature_name = gene, plotTitle = paste("Native_1 -", gene))
  p2 <- plotRaster(rastList$Native_2, feature_name = gene, plotTitle = paste("Native_2 -", gene))
  ggsave(file.path(cfg$outdir, paste0(gene, "_raster.png")), p1 + p2, width = 10, height = 5, dpi = 300)

  # Correlation plots
  if (gene %in% rownames(sc)) {
    p_cor <- plotCorrelationGeneExp(
      list(Native_1 = rastList$Native_1, Native_2 = rastList$Native_2), sc, gene
    )
    ggsave(file.path(cfg$outdir, paste0(gene, "_correlation.png")), p_cor, width = 10, height = 5, dpi = 300)
  }

  # Linear regression
  p_lr <- linearRegression(input = ss, gene = gene)
  ggsave(file.path(cfg$outdir, paste0(gene, "_linearRegression.png")), p_lr, width = 10, height = 5, dpi = 300)

  # Pixel class
  p_pc <- pixelClass(input = ss, gene = gene)
  ggsave(file.path(cfg$outdir, paste0(gene, "_pixelClass.png")), p_pc, width = 10, height = 5, dpi = 300)

  cat("  Saved plots for:", gene, "\n")
}

cat("\n=== Done. All outputs written to:", cfg$outdir, "===\n")