#!/usr/bin/env Rscript
# aggregates the per-pair Overall_Similarity.csv files produced by STcompare.R into a matrix heatmap

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(ggplot2)
  library(reshape2)
})

p <- arg_parser("Print pairwise STcompare similarity results as summary")
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

# reading Overall_Similarity.csv
read_one <- function(dir) {
  f <- file.path(dir, "Results", "Overall_Similarity.csv")
  if (!file.exists(f)) {
    warning("Missing Overall_Similarity.csv in ", dir, " -- skipping")
    return(NULL)
  }
  read.csv(f, stringsAsFactors = FALSE)
}

pair_results <- do.call(rbind, lapply(pair_dirs, read_one))
if (is.null(pair_results) || nrow(pair_results) == 0) {
  stop("No Overall_Similarity.csv files could be read. Rrerun compare first.")
}

# 1. line-based similarity table (CSV)
write.csv(pair_results, file.path(argv$outdir, "All_Pairwise_Similarity_Long.csv"), row.names = FALSE)

# 2. sample matrix (CSV)
#   Samples are ordered by their DonorX_dayY_difZ_orgW components 
#   so that heatmap lists organoids in the same order as the 
#   reference similarity heatmap
#   grouped by day first (120, then 70, then 40), then by donor,
#   then dif, then org
day_levels   <- c("day120", "day70", "day40")
donor_levels <- c("Donor2", "Donor1")
dif_levels   <- c("dif2", "dif1")
org_levels   <- c("org2", "org1")

order_samples <- function(sample_names) {
  m <- regmatches(
    sample_names,
    regexec("^(Donor[0-9]+)_(day[0-9]+)_(dif[0-9]+)_(org[0-9]+)$", sample_names)
  )
  parsed <- do.call(rbind, lapply(seq_along(sample_names), function(i) {
    parts <- m[[i]]
    if (length(parts) != 5) {
      stop("Sample name does not match expected DonorX_dayY_difZ_orgW pattern: ", sample_names[i])
    }
    data.frame(
      sample = sample_names[i], donor = parts[2], day = parts[3],
      dif = parts[4], org = parts[5], stringsAsFactors = FALSE
    )
  }))
  parsed$donor <- factor(parsed$donor, levels = donor_levels)
  parsed$day   <- factor(parsed$day,   levels = day_levels)
  parsed$dif   <- factor(parsed$dif,   levels = dif_levels)
  parsed$org   <- factor(parsed$org,   levels = org_levels)
  parsed[order(parsed$day, parsed$donor, parsed$dif, parsed$org), ]$sample
}

# top-to-bottom display order, matching the reference heatmap, used directly for the CSV matrix 
samples <- order_samples(unique(c(pair_results$sample_aligned, pair_results$sample_reference)))
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

# 3. heatmap (PDF)
plot_levels <- rev(samples)
mat_df <- melt(mat, varnames = c("sample1", "sample2"), value.name = "similarity")
mat_df$sample1 <- factor(mat_df$sample1, levels = plot_levels)
mat_df$sample2 <- factor(mat_df$sample2, levels = plot_levels)
mat_df <- mat_df[as.integer(mat_df$sample1) >= as.integer(mat_df$sample2), ]

# colour scale is capped to the highest real similarity value outside of the diagonal
offdiag_vals <- mat[row(mat) != col(mat)]
max_val <- max(offdiag_vals, na.rm = TRUE)

# label colour flips to black on the brightest (most yellow) tiles so txt is readable
label_threshold <- 0.7 * max_val
mat_df$label_color <- ifelse(is.na(mat_df$similarity), NA,
  ifelse(mat_df$similarity > label_threshold, "black", "white"))

heatmap_plot <- ggplot(mat_df, aes(x = sample1, y = sample2, fill = similarity)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(similarity), "", sprintf("%.2f", similarity)), color = label_color), size = 3) +
  scale_color_identity() +
  scale_fill_viridis_c(option = "magma", limits = c(0, max_val), oob = scales::squish, name = "Similarity", na.value = "white") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.3),
    axis.title = element_blank()
  )

ggsave(
  file.path(argv$outdir, "spatial_similarity_heatmap.pdf"),
  heatmap_plot, width = 10, height = 9, dpi = 300, bg = "white"
)

print(paste0("Aggregated", nrow(pair_results), "pairs across", length(samples), "samples.\n"))