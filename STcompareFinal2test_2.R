# Library loading
library(STcompare)
library(SpatialExperiment)
library(SEraster)
library(hdf5r)
library(Seurat) 
library(patchwork)
library(ggplot2)
library(jsonlite)

# Data loading
path_counts1 <- "/Users/adrhovska/Desktop/STdata/Native_1_ST/filtered_feature_bc_matrix.h5"
path_counts2 <- "/Users/adrhovska/Desktop/STdata/Native_2_ST/filtered_feature_bc_matrix.h5"
path_pos1_aligned <- "/Users/adrhovska/Desktop/STdata/Aligned_Native_1_2/Native_1_aligned_to_Native_2_WITH_BARCODES.csv"
path_pos2_spatial_dir <- "/Users/adrhovska/Desktop/STdata/Native_2_ST/spatial"

# Check for presence of Gene expression object 
get_gene_expression <- function(counts) {
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  counts
}

# Read the files with assigned barcodes (check prior) - aligned ones 
read_pos_with_barcodes <- function(path) {
  pos <- read.csv(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Check for required columns (continues if FALSE)
  if (!all(c("barcode", "x", "y") %in% colnames(pos))) {
    stop("Aligned file does not contain columns: barcode, x, y")
  }
  
  coord <- data.frame(
    x = as.numeric(pos$x),
    y = as.numeric(pos$y),
    row.names = pos$barcode
  )
  
  if (any(!is.finite(coord$x)) || any(!is.finite(coord$y))) {
    stop("coordinates contain non-finite values.")
  }
  
  coord
}

# Read positions which are not aligned and construct paths
read_pos <- function(spatial_dir, scale_type = "hires") {
  
  pos_path <- file.path(spatial_dir, "tissue_positions.csv")
  scale_path <- file.path(spatial_dir, "scalefactors_json.json")
  
  pos <- read.csv(
    pos_path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!all(c("pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(pos))) {
    stop("tissue_positions.csv does not contain the expected columns") # Check 
  }
  
  # Read and choose scale factor 
  scales <- jsonlite::fromJSON(scale_path)
  
  if (scale_type == "hires") {
    scale_factor <- scales$tissue_hires_scalef
  } else if (scale_type == "lowres") {
    scale_factor <- scales$tissue_lowres_scalef
  } else {
    stop("scale_type must be either 'hires' or 'lowres'.")
  }
  
cat(
  "Using", scale_type,
  "scale factor:", scale_factor, "\n"
)
# Scaling 
  coords <- data.frame(
    x = as.numeric(pos$pxl_col_in_fullres) * scale_factor,
    y = as.numeric(pos$pxl_row_in_fullres) * scale_factor,
    row.names = rownames(pos)
  )
  
  if (any(!is.finite(coords$x)) || any(!is.finite(coords$y))) {
    stop("scaled coordinates contain non-finite values") # Check 
  }
  
  coords
}
# Match gene expression counts to spatial coordinate matrix
match_counts_to_positions <- function(counts, pos, sample_name) {
  exact_common <- intersect(colnames(counts), rownames(pos))
  cat(sample_name, "barcode matches:", length(exact_common), "\n")
  if (length(exact_common) > 0) {
    counts_matched <- counts[, exact_common, drop = FALSE]
    pos_matched <- pos[exact_common, , drop = FALSE]
  } else {
    # cleaning 
    count_key <- sub("-1$", "", colnames(counts))
    pos_key <- sub("-1$", "", rownames(pos))
    common_key <- intersect(count_key, pos_key)
    cat(sample_name, "cleaned barcode matches:", length(common_key), "\n")
    if (length(common_key) == 0) {
      stop(paste0("No matching barcodes found for ", sample_name))
    }
    
    count_idx <- match(common_key, count_key)
    pos_idx <- match(common_key, pos_key)
    
    counts_matched <- counts[, count_idx, drop = FALSE]
    pos_matched <- pos[pos_idx, , drop = FALSE]
    
    rownames(pos_matched) <- colnames(counts_matched)
  }
  
  coords <- cbind(
    x = as.numeric(pos_matched$x),
    y = as.numeric(pos_matched$y)
  )
  
  rownames(coords) <- rownames(pos_matched)
  
  if (ncol(counts_matched) == 0 || nrow(coords) == 0) {
    stop(paste0(sample_name, " has zero matched spots."))
  }
  
  if (!all(is.finite(coords))) {
    stop(paste0(sample_name, " has non-finite coordinates."))
  }
  if (!identical(colnames(counts_matched), rownames(coords))) {
  stop(paste0(sample_name, " count columns and coordinate rows are not aligned."))
  }
  
  cat(sample_name, "final matched spots:", ncol(counts_matched), "\n")
  
  list(
    counts = counts_matched,
    coords = coords
  )
}
# Range check for coordinates 
coordinate_ranges_check <- function(coords1, coords2) {
  coord_ranges <- data.frame(
    sample = c("Native_1_aligned", "Native_2_scaled"),
    min_x = c(min(coords1[, "x"]), min(coords2[, "x"])),
    max_x = c(max(coords1[, "x"]), max(coords2[, "x"])),
    min_y = c(min(coords1[, "y"]), min(coords2[, "y"])),
    max_y = c(max(coords1[, "y"]), max(coords2[, "y"]))
  )
  
  print(coord_ranges)
  
  x_overlap <- max(coords1[, "x"]) >= min(coords2[, "x"]) &&
    max(coords2[, "x"]) >= min(coords1[, "x"])
  
  y_overlap <- max(coords1[, "y"]) >= min(coords2[, "y"]) &&
    max(coords2[, "y"]) >= min(coords1[, "y"])
  
  if (!x_overlap || !y_overlap) {
    stop("No overlap")
  }
}

# Data read 
counts1 <- get_gene_expression(Read10X_h5(path_counts1))
counts2 <- get_gene_expression(Read10X_h5(path_counts2))

cat("\nNative 1 count matrix dimensions:\n") # Have to replace the Natives to the actual sample names in the future
print(dim(counts1))

cat("\nNative 2 count matrix dimensions:\n")
print(dim(counts2))


pos1_scaled <- read_pos_with_barcodes(path_pos1_aligned)

pos2_scaled <- read_pos(
  spatial_dir = path_pos2_spatial_dir,
  scale_type = "hires"
)

matched1 <- match_counts_to_positions(counts1, pos1_scaled, "Native_1")
matched2 <- match_counts_to_positions(counts2, pos2_scaled, "Native_2")

counts1_matched <- matched1$counts
coords1 <- matched1$coords

counts2_matched <- matched2$counts
coords2 <- matched2$coords


epithelial_genes      <- c("KRT4", "KRT5", "IVL")
smooth_muscle_genes   <- c("SMTN", "CALD1", "CSRP1", "TAGLN")
skeletal_muscle_genes <- c("TNNC1", "TNNC2", "ACTC1", "MYH8")

genes_of_interest <- c(
  epithelial_genes,
  smooth_muscle_genes,
  skeletal_muscle_genes
)

genes_of_interest <- genes_of_interest[
  genes_of_interest %in% rownames(counts1_matched) &
    genes_of_interest %in% rownames(counts2_matched)
]

if (length(genes_of_interest) == 0) {
  stop("None of the genes of interest are present in both count matrices.")
}

cat("\nGenes used:\n")
print(genes_of_interest)

counts1_matched <- counts1_matched[genes_of_interest, , drop = FALSE]
counts2_matched <- counts2_matched[genes_of_interest, , drop = FALSE]


spe1 <- SpatialExperiment(
  assays = list(counts = counts1_matched),
  spatialCoords = coords1
)

spe2 <- SpatialExperiment(
  assays = list(counts = counts2_matched),
  spatialCoords = coords2
)

df_coords <- rbind(
  data.frame(
    x = coords1[, "x"],
    y = coords1[, "y"],
    sample = "Native_1_aligned"
  ),
  data.frame(
    x = coords2[, "x"],
    y = coords2[, "y"],
    sample = "Native_2_scaled"
  )
)

print(
  ggplot(df_coords, aes(x = x, y = y, colour = sample)) +
    geom_point(size = 0.4, alpha = 0.5) +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Coordinate overlap check")
)

rastList <- SEraster::rasterizeGeneExpression(
  list(Native_1 = spe1, Native_2 = spe2),
  assay_name = "counts",
  resolution = 150,
  square = FALSE
)

all_x <- c(coords1[, "x"], coords2[, "x"])
all_y <- c(coords1[, "y"], coords2[, "y"])

x_pad <- diff(range(all_x)) * 0.05
y_pad <- diff(range(all_y)) * 0.05

shared_xlim <- range(all_x) + c(-x_pad, x_pad)
shared_ylim <- range(all_y) + c(-y_pad, y_pad)

sc <- spatialCorrelationGeneExp(
  list(Native_1 = rastList$Native_1, Native_2 = rastList$Native_2),
  nThreads = 4
)

ss <- spatialSimilarity(
  list(Native_1 = rastList$Native_1, Native_2 = rastList$Native_2)
)

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(sc)]

results <- sc[
  genes_of_interest,
  c("correlationCoef", "pValuePermuteX", "pValuePermuteY"),
  drop = FALSE
]

results$empirical_pval <- pmax(
  results$pValuePermuteX,
  results$pValuePermuteY
)

results$cell_type <- ifelse(
  rownames(results) %in% epithelial_genes,
  "Epithelial",
  ifelse(
    rownames(results) %in% smooth_muscle_genes,
    "Smooth Muscle",
    "Skeletal Muscle"
  )
)

print(results)

# Plots 
for (gene in genes_of_interest) {
  p1 <- plotRaster(
    rastList$Native_1,
    feature_name = gene,
    plotTitle = paste("Native_1 aligned -", gene)
  )
  
  p2 <- plotRaster(
    rastList$Native_2,
    feature_name = gene,
    plotTitle = paste("Native_2 scaled -", gene)
  )
  
  print(p1 + p2)
}


for (gene in genes_of_interest) {
  print(plotCorrelationGeneExp(
    list(Native_1 = rastList$Native_1, Native_2 = rastList$Native_2),
    sc,
    gene
  ))
}


plot_dir <- "/Users/adrhovska/Desktop/STdata/STcompare_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


for (gene in genes_of_interest) {
  
  p1 <- plotRaster(
    rastList$Native_1,
    feature_name = gene,
    plotTitle = paste("Native_1 aligned -", gene)
  )
  
  p2 <- plotRaster(
    rastList$Native_2,
    feature_name = gene,
    plotTitle = paste("Native_2 scaled -", gene)
  )
  
  p <- p1 + p2
  
  ggsave(
    filename = file.path(plot_dir, paste0(gene, "_raster.png")),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300
  )
}

for (gene in genes_of_interest) {
  
  p <- plotCorrelationGeneExp(
    list(Native_1 = rastList$Native_1, Native_2 = rastList$Native_2),
    sc,
    gene
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0(gene, "_correlation.png")),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300
  )
}

for (gene in genes_of_interest) {
  
  p_lr <- linearRegression(input = ss, gene = gene)
  
  ggsave(
    filename = file.path(plot_dir, paste0(gene, "_linearRegression.png")),
    plot = p_lr,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  p_pc <- pixelClass(input = ss, gene = gene)
  
  ggsave(
    filename = file.path(plot_dir, paste0(gene, "_pixelClass.png")),
    plot = p_pc,
    width = 10,
    height = 5,
    dpi = 300
  )
}