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

# Extracting from command line arguments and parsing them
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    counts1 = NULL,
    counts2 = NULL,
    pos1    = NULL,
    spatial2 = NULL,
    outdir  = "./STcompare_out",
    scale   = "hires",
    res     = 150,
    threads = 4,
    sample_aligned = "Sample_1",
    sample_reference = "Sample_2"
  )

  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    val <- if (i + 1 <= length(args)) args[i + 1] else stop(paste0("Missing value for --", key))
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
    cat("    --counts1 \\\n")
    cat("    --counts2 \\\n")
    cat("    --pos1 \\\n")
    cat("    --spatial2 \\\n")
    cat("    [--outdir  ./STcompare_out] \\\n") # Those in [] mean that they are preset as a default value
    cat("    [--scale   hires|lowres] \\\n")
    cat("    [--res     150] \\\n")
    cat("    [--threads 4]\n")
    quit(status = 1)
  }

  defaults
}

cfg <- parse_args(args)

sample_aligned_name <- cfg$sample_aligned
sample_reference_name <- cfg$sample_reference

cat("counts1  :", cfg$counts1,  "\n")
cat("counts2  :", cfg$counts2,  "\n")
cat("pos1     :", cfg$pos1,     "\n")
cat("spatial2 :", cfg$spatial2, "\n")
cat("outdir   :", cfg$outdir,   "\n")
cat("scale    :", cfg$scale,    "\n")
cat("res      :", cfg$res,      "\n")
cat("threads  :", cfg$threads,  "\n\n")
cat("sample_aligned  :", cfg$sample_aligned,  "\n")
cat("sample_reference  :", cfg$sample_reference,  "\n")

dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)

comparison_name <- paste0(
  sample_aligned_name,
  "_vs_",
  sample_reference_name,
  "_",
  cfg$scale,
  "_res",
  cfg$res
)

dir_comparison <- file.path(cfg$outdir, comparison_name)

dir.create(dir_comparison, showWarnings = FALSE, recursive = TRUE)

dir_results <- file.path(dir_comparison, "00_results")
dir_qc <- file.path(dir_comparison, "01_coordinate_qc")
dir_raster <- file.path(dir_comparison, "02_raster_plots")
dir_correlation <- file.path(dir_comparison, "03_correlation_plots")
dir_linear <- file.path(dir_comparison, "04_linear_regression")
dir_pixel <- file.path(dir_comparison, "05_pixel_class")

output_dirs <- c(
  dir_results,
  dir_qc,
  dir_raster,
  dir_correlation,
  dir_linear,
  dir_pixel
)

invisible(lapply(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

epithelial_genes      <- c("KRT4", "KRT5", "IVL")
smooth_muscle_genes   <- c("SMTN", "CALD1", "CSRP1", "TAGLN")
skeletal_muscle_genes <- c("TNNC1", "TNNC2", "ACTC1", "MYH8")

genes_of_interest <- c(epithelial_genes, smooth_muscle_genes, skeletal_muscle_genes)

get_gene_expression <- function(counts) {
  if (!is.list(counts)) {
    return(counts)
  }
  if ("Gene Expression" %in% names(counts)) {
    return(counts[["Gene Expression"]])
  }
  return (counts[[1]])
}

read_aligned_positions <- function(path, sample_name = "sample") {
  
  pos <- read.csv(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  required_cols <- c("barcode", "x", "y")
  
  if (!all(required_cols %in% colnames(pos))) {
    stop(paste0(
      sample_name,
      " aligned file missing expected columns: ",
      paste(required_cols, collapse = ", ")
    ))
  }
  
  coord <- data.frame(
    x = as.numeric(pos$x),
    y = as.numeric(pos$y),
    row.names = pos$barcode
  )
  
  if (any(!is.finite(coord$x)) || any(!is.finite(coord$y))) {
    stop(paste0(sample_name, " aligned coordinates contain non-finite values."))
  }
  
  coord
}

read_visium_positions <- function(
  spatial_dir,
  scale_type = "hires",
  sample_name = "sample"
) {
  
  pos_path <- file.path(spatial_dir, "tissue_positions.csv")
  scale_path <- file.path(spatial_dir, "scalefactors_json.json")
  
  pos <- read.csv(
    pos_path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  required_cols <- c("pxl_row_in_fullres", "pxl_col_in_fullres")
  
  if (!all(required_cols %in% colnames(pos))) {
    stop(paste0(
      sample_name,
      " tissue_positions.csv missing expected columns: ",
      paste(required_cols, collapse = ", ")
    ))
  }
  
  scales <- jsonlite::fromJSON(scale_path)
  
  if (scale_type == "hires") {
    scale_factor <- scales$tissue_hires_scalef
  } else if (scale_type == "lowres") {
    scale_factor <- scales$tissue_lowres_scalef
  } else {
    stop("scale_type must be either 'hires' or 'lowres'.")
  }
  
  cat(sample_name, scale_type, "scale factor:", scale_factor, "\n")
  
  pos_scaled <- data.frame(
    x = as.numeric(pos$pxl_col_in_fullres) * scale_factor,
    y = as.numeric(pos$pxl_row_in_fullres) * scale_factor,
    row.names = rownames(pos)
  )
  
  if (any(!is.finite(pos_scaled$x)) || any(!is.finite(pos_scaled$y))) {
    stop(paste0(sample_name, " scaled coordinates contain non-finite values."))
  }
  
  pos_scaled
}

# Coord check
check_coordinate_system <- function(coords1, coords2, sample_aligned_name, sample_reference_name) {
  
  ranges <- data.frame(
    sample = c(sample_aligned_name, sample_reference_name),
    min_x = c(min(coords1[, "x"]), min(coords2[, "x"])),
    max_x = c(max(coords1[, "x"]), max(coords2[, "x"])),
    min_y = c(min(coords1[, "y"]), min(coords2[, "y"])),
    max_y = c(max(coords1[, "y"]), max(coords2[, "y"]))
  )
  
  ranges$width <- ranges$max_x - ranges$min_x
  ranges$height <- ranges$max_y - ranges$min_y
  
  print(ranges)
  
  width_ratio <- max(ranges$width) / min(ranges$width)
  height_ratio <- max(ranges$height) / min(ranges$height)
  
  cat("Width ratio:", width_ratio, "\n")
  cat("Height ratio:", height_ratio, "\n")
  
  x_overlap <- max(coords1[, "x"]) >= min(coords2[, "x"]) &&
    max(coords2[, "x"]) >= min(coords1[, "x"])
  
  y_overlap <- max(coords1[, "y"]) >= min(coords2[, "y"]) &&
    max(coords2[, "y"]) >= min(coords1[, "y"])
  
  if (!x_overlap || !y_overlap) {
    warning(
      "Coordinate ranges do not overlap, samples are not in the same aligned coordinate space."
    )
  }
  
  if (width_ratio > 3 || height_ratio > 3) {
    warning(
      "Coordinate ranges differ strongly in scale, one sample may be fullres while the other is hires/lowres."
    )
  }
  
  invisible(NULL)
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

counts1 <- get_gene_expression(Read10X_h5(cfg$counts1))
counts2 <- get_gene_expression(Read10X_h5(cfg$counts2))
cat(sample_aligned_name, "dims:", dim(counts1), "\n")
cat(sample_reference_name, "dims:", dim(counts2), "\n")

pos1 <- read_aligned_positions(
  path = cfg$pos1,
  sample_name = sample_aligned_name
)

pos2 <- read_visium_positions(
  spatial_dir = cfg$spatial2,
  scale_type = cfg$scale,
  sample_name = sample_reference_name
)

check_coordinate_system(
  coords1 = pos1,
  coords2 = pos2,
  sample_aligned_name = sample_aligned_name,
  sample_reference_name = sample_reference_name
)

matched1 <- match_counts_to_positions(counts1, pos1, sample_aligned_name)
matched2 <- match_counts_to_positions(counts2, pos2, sample_reference_name)

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
  data.frame(
    x = coords1[, "x"],
    y = coords1[, "y"],
    sample = sample_aligned_name
  ),
  data.frame(
    x = coords2[, "x"],
    y = coords2[, "y"],
    sample = sample_reference_name
  )
)

p_overlap <- ggplot(df_coords, aes(x = x, y = y, colour = sample)) +
  geom_point(size = 0.4, alpha = 0.5) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Coordinate overlap check",
    x = "x coordinate",
    y = "y coordinate",
    colour = "Sample"
  )

ggsave(
  file.path(dir_qc, "Coordinate_Overlap.png"),
  p_overlap,
  width = 7,
  height = 6,
  dpi = 200
)

# Object building 
spe1 <- SpatialExperiment(
  assays = list(counts = counts1_matched),
  spatialCoords = coords1
)

spe2 <- SpatialExperiment(
  assays = list(counts = counts2_matched),
  spatialCoords = coords2
)

spe_list <- setNames(
  list(spe1, spe2),
  c(sample_aligned_name, sample_reference_name)
)

# Rasterization
rastList <- SEraster::rasterizeGeneExpression(
  spe_list,
  assay_name = "counts",
  resolution = cfg$res,
  square = FALSE
)

# STcompare
sc <- spatialCorrelationGeneExp(
  rastList,
  nThreads = cfg$threads
)

ss <- spatialSimilarity(
  rastList
)

genes_in_sc <- genes_of_interest[
  genes_of_interest %in% rownames(sc)
]

results <- sc[
  genes_in_sc,
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

write.csv(
  results,
  file.path(dir_results, "Results_Table.csv")
)

# Identify raster assay name
assays1 <- SummarizedExperiment::assayNames(rastList[[sample_aligned_name]])
assays2 <- SummarizedExperiment::assayNames(rastList[[sample_reference_name]])

common_assays <- intersect(assays1, assays2)

if (length(common_assays) == 0) {
  stop("No shared assay names between rasterised samples.")
}

rast_assay <- if ("counts" %in% common_assays) {
  "counts"
} else {
  common_assays[1]
}

cat("Using raster assay:", rast_assay, "\n")


# Shared spatial limits
all_x <- c(coords1[, "x"], coords2[, "x"])
all_y <- c(coords1[, "y"], coords2[, "y"])

all_x <- all_x[is.finite(all_x)]
all_y <- all_y[is.finite(all_y)]

raster_resolution <- cfg$res
raster_padding_multiplier <- 2

spatial_pad <- raster_resolution * raster_padding_multiplier

shared_xlim <- c(min(all_x) - spatial_pad, max(all_x) + spatial_pad)
shared_ylim <- c(min(all_y) - spatial_pad, max(all_y) + spatial_pad)

spatial_unit_label <- paste0(cfg$scale, " image pixels")
expression_unit_label <- "rasterised raw counts"


get_shared_gene_limits <- function(rastList, gene, assay_name, name1, name2) {
  
  vals1 <- as.numeric(
    SummarizedExperiment::assay(
      rastList[[name1]],
      assay_name
    )[gene, ]
  )
  
  vals2 <- as.numeric(
    SummarizedExperiment::assay(
      rastList[[name2]],
      assay_name
    )[gene, ]
  )
  
  vals <- c(vals1, vals2)
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) {
    stop(paste0("No finite raster values found for gene ", gene))
  }
  
  limits <- range(vals, na.rm = TRUE)
  
  if (limits[1] == limits[2]) {
    limits <- limits + c(-0.5, 0.5)
  }
  
  limits
}

make_single_raster <- function(rast, name, gene, gene_limits,
                                rast_assay, shared_xlim, shared_ylim,
                                spatial_unit_label,
                                expression_unit_label) {
  SEraster::plotRaster(
    rast,
    assay_name  = rast_assay,
    feature_name = gene,
    plotTitle   = paste(name, "-", gene)
  ) +
    ggplot2::scale_fill_viridis_c(
      limits = gene_limits,
      oob    = scales::squish,
      name   = paste0(gene, "\n", expression_unit_label)
    ) +
    ggplot2::coord_sf(
      xlim   = shared_xlim,
      ylim   = shared_ylim,
      expand = FALSE,
      clip   = "off"
    ) +
    ggplot2::labs(
      x = paste0("x coordinate (", spatial_unit_label, ")"),
      y = paste0("y coordinate (", spatial_unit_label, ")")
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10))
}

make_raster_pair <- function(gene, rastList, rast_assay,
                              name1, name2,
                              shared_xlim, shared_ylim,
                              spatial_unit_label,
                              expression_unit_label) {
  
  gene_limits <- get_shared_gene_limits(
    rastList = rastList,
    gene = gene,
    assay_name = rast_assay,
    name1 = name1,
    name2 = name2
  )
  
  p1 <- make_single_raster(
    rast                 = rastList[[name1]],
    name                 = name1,
    gene                 = gene,
    gene_limits          = gene_limits,
    rast_assay           = rast_assay,
    shared_xlim          = shared_xlim,
    shared_ylim          = shared_ylim,
    spatial_unit_label   = spatial_unit_label,
    expression_unit_label = expression_unit_label
)

p2 <- make_single_raster(
    rast                 = rastList[[name2]],
    name                 = name2,
    gene                 = gene,
    gene_limits          = gene_limits,
    rast_assay           = rast_assay,
    shared_xlim          = shared_xlim,
    shared_ylim          = shared_ylim,
    spatial_unit_label   = spatial_unit_label,
    expression_unit_label = expression_unit_label
)
  
  p1 + p2 +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
}

for (gene in genes_of_interest) {
  
  if (gene %in% genes_in_sc) {
    
    p_raster <- make_raster_pair(
      gene                 = gene,
      rastList             = rastList,
      rast_assay           = rast_assay,
      name1                = sample_aligned_name,
      name2                = sample_reference_name,
      shared_xlim          = shared_xlim,
      shared_ylim          = shared_ylim,
     spatial_unit_label   = spatial_unit_label,
      expression_unit_label = expression_unit_label
)
    
    ggsave(
      file.path(dir_raster, paste0(gene, "_Raster.png")),
      p_raster,
      width = 11,
      height = 5.5,
      dpi = 300,
      bg = "white"
    )
    
    p_cor <- plotCorrelationGeneExp(
      rastList,
      sc,
      gene
    )
    
    ggsave(
      file.path(dir_correlation, paste0(gene, "_Correlation.png")),
      p_cor,
      width = 10,
      height = 5,
      dpi = 300,
      bg = "white"
    )
  }
  
  p_lr <- linearRegression(input = ss, gene = gene)
  
  ggsave(
    file.path(dir_linear, paste0(gene, "_LinearRegression.png")),
    p_lr,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  p_pc <- pixelClass(input = ss, gene = gene)
  
  ggsave(
    file.path(dir_pixel, paste0(gene, "_PixelClass.png")),
    p_pc,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"
  )

}