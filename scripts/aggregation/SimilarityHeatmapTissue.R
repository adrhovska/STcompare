# build_tissue_cluster_heatmaps.R
#
# Aggregates the Cluster_Overall_Similarity.csv files produced by
# STcompare_with_tissue_clusters.R into one
# all-organoids similarity heatmap per tissue-type cluster 

suppressPackageStartupMessages({
  library(ggplot2)
  library(argparser, quietly = TRUE)
})

p <- arg_parser("build_tissue_cluster_heatmaps")
p <- add_argument(p, "--stcompare_root", help = "Directory containing one subdirectory per STcompare sample-pair run")
p <- add_argument(p, "--outdir", help = "Where to save the heatmap PDFs", default = "./Cluster_Heatmaps")
argv <- parse_args(p)

if (is.na(argv$stcompare_root)) {
  print("ERROR: --stcompare_root is required")
  quit(status = 1)
}
dir.create(argv$outdir, showWarnings = FALSE, recursive = TRUE)

# collecting every Cluster_Overall_Similarity.csv under the root directory
csv_paths <- list.files(argv$stcompare_root,
  pattern = "^Cluster_Overall_Similarity\\.csv$",
  recursive = TRUE, full.names = TRUE
)
if (length(csv_paths) == 0) {
  print(paste("ERROR: No Cluster_Overall_Similarity.csv files found under", argv$stcompare_root))
  quit(status = 1)
}
cat("Found", length(csv_paths), "pairwise result file(s)\n")

all_pairs <- do.call(rbind, lapply(csv_paths, read.csv, stringsAsFactors = FALSE))

# fixing ordering
all_samples <- sort(unique(c(all_pairs$sample_aligned, all_pairs$sample_reference)))

# building a symmetric similarity matrix for one cluster
#   @df: subset of all_pairs for a single cluster
#   @samples: full ordered vector of sample names (rows/cols of the matrix)
make_similarity_matrix <- function(df, samples) {
  m <- matrix(NA_real_, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
  diag(m) <- 1
  for (i in seq_len(nrow(df))) {
    a <- df$sample_aligned[i]
    b <- df$sample_reference[i]
    if (a %in% samples && b %in% samples) {
      val <- df$mean_percent_similarity[i] / 100
      m[a, b] <- val
      m[b, a] <- val
    }
  }
  m
}

  samples <- rownames(m)
  df <- expand.grid(
    sample_row = factor(samples, levels = samples),
    sample_col = factor(samples, levels = samples)
  )
  df$value <- as.vector(m[cbind(as.character(df$sample_row), as.character(df$sample_col))])
  row_idx <- match(df$sample_row, samples)
  col_idx <- match(df$sample_col, samples)
  df[row_idx <= col_idx, ]
}

plot_similarity_heatmap <- function(df_long, title) {
  ggplot(df_long, aes(x = sample_col, y = sample_row, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value))), size = 2.6) +
    scale_fill_viridis_c(option = "inferno", limits = c(0, 1), name = "Similarity", na.value = "grey90") +
    scale_x_discrete(position = "bottom") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle(title)
}

clusters <- unique(all_pairs$cluster)
for (cl in clusters) {
  df_cl <- all_pairs[all_pairs$cluster == cl, ]
  m <- make_similarity_matrix(df_cl, all_samples)
  df_long <- matrix_to_upper_tri_long(m)
  plt <- plot_similarity_heatmap(df_long, paste("Spatial similarity -", cl))
  out_path <- file.path(argv$outdir, paste0(cl, "_spatial_similarity_heatmap.pdf"))
  ggsave(out_path, plt, width = 10, height = 8, dpi = 300, bg = "white")
  cat("Saved", out_path, "\n")
}

print(paste0(Done.))