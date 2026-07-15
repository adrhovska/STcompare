#!/usr/bin/env Rscript
# aggregates the per-pair Overall_Similarity.csv files produced by STcompare.R

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(ggplot2)
  library(reshape2)
})

p <- arg_parser("Print pairwise STcompare similarity results as one summary")
p <- add_argument(p, "--project_dir", help = "Project directory containing the pairwise output folder", default = ".")
p <- add_argument(p, "--pairwise_dir", help = "Name of the pairwise output directory", default = "STcompare_pairwise")
p <- add_argument(p, "--outdir", help = "Where to write the aggregated outputs", default = "./STcompare_summary")
argv <- parse_args(p)

dir.create(argv$outdir, showWarnings = FALSE, recursive = TRUE)

pairwise_root <- file.path(argv$project_dir, argv$pairwise_dir)
pair_dirs <- list.dirs(pairwise_root, recursive = FALSE)

if (length(pair_dirs) == 0) {
  stop("No pairwise output directories found in ", pairwise_root)
}

# reading one pair's Overall_Similarity.csv
read_one <- function(dir) {
  f <- file.path(dir, "Results", "Overall_Similarity.csv")
  if (!file.exists(f)) {
    warning("Missing Overall_Similarity.csv in ", dir, " -- skipping (rerun STcompare.R with the snippet added)")
    return(NULL)
  }
  read.csv(f, stringsAsFactors = FALSE)
}

pair_results <- do.call(rbind, lapply(pair_dirs, read_one))
if (is.null(pair_results) || nrow(pair_results) == 0) {
  stop("No Overall_Similarity.csv files could be read. Add the snippet to STcompare.R and rerun star_compare first.")
}

# 1. long-format table: one row per unique pair
write.csv(pair_results, file.path(argv$outdir, "All_Pairwise_Similarity_Long.csv"), row.names = FALSE)

# 2. full symmetric sample x sample matrix, diagonal = 1
samples <- sort(unique(c(pair_results$sample_aligned, pair_results$sample_reference)))
mat <- matrix(NA_real_, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
diag(mat) <- 1

for (i in seq_len(nrow(pair_results))) {
  a <- pair_results$sample_aligned[i]
  b <- pair_results$sample_reference[i]
  v <- pair_results$mean_percent_similarity[i]
  mat[a, b] <- v
  mat[b, a] <- v
}

write.csv(mat, file.path(argv$outdir, "All_Pairwise_Similarity_Matrix.csv"), row.names = TRUE)

# 3. triangular heatmap (upper triangle, values printed in each cell)
mat_df <- melt(mat, varnames = c("sample1", "sample2"), value.name = "similarity")
mat_df$sample1 <- factor(mat_df$sample1, levels = samples)
mat_df$sample2 <- factor(mat_df$sample2, levels = samples)
mat_df <- mat_df[as.integer(mat_df$sample1) <= as.integer(mat_df$sample2), ]

heatmap_plot <- ggplot(mat_df, aes(x = sample1, y = sample2, fill = similarity)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(similarity), "", sprintf("%.2f", similarity))), size = 3) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), name = "Similarity", na.value = "white") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  file.path(argv$outdir, "spatial_similarity_heatmap.pdf"),
  heatmap_plot, width = 10, height = 9, dpi = 300, bg = "white"
)

cat("Aggregated", nrow(pair_results), "pairs across", length(samples), "samples.\n")