
aligned_path <- "/Users/adrhovska/Desktop/STdata/Aligned_Native_1_2/Native_1_aligned_to_Native_2.csv"

native1_pos_path <- "/Users/adrhovska/Desktop/STdata/Native_1_ST/spatial/tissue_positions.csv"

output_path <- "/Users/adrhovska/Desktop/STdata/Aligned_Native_1_2/Native_1_aligned_to_Native_2_WITH_BARCODES.csv"


aligned <- read.csv(
  aligned_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Aligned file columns:\n")
print(colnames(aligned))
cat("Aligned rows:", nrow(aligned), "\n")

# If the file has an accidental index column like X, ignore it.
if (!all(c("x", "y") %in% colnames(aligned))) {
  stop("Aligned file does not contain columns called x and y.")
}

aligned_xy <- data.frame(
  x = as.numeric(aligned$x),
  y = as.numeric(aligned$y)
)

aligned_xy <- aligned_xy[
  is.finite(aligned_xy$x) & is.finite(aligned_xy$y),
  ,
  drop = FALSE
]

cat("Aligned finite coordinate rows:", nrow(aligned_xy), "\n")


native1_pos <- read.csv(
  native1_pos_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Handle old 10x format with no header
if (!"barcode" %in% colnames(native1_pos)) {
  native1_pos <- read.csv(
    native1_pos_path,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(native1_pos) <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )
}

cat("Native 1 tissue_positions columns:\n")
print(colnames(native1_pos))

native1_pos <- native1_pos[native1_pos$in_tissue == 1, , drop = FALSE]

cat("Native 1 in-tissue rows:", nrow(native1_pos), "\n")



if (nrow(aligned_xy) != nrow(native1_pos)) {
  stop(
    paste0(
      "Cannot safely attach barcodes.\n",
      "Aligned coordinate rows: ", nrow(aligned_xy), "\n",
      "Native 1 in-tissue rows: ", nrow(native1_pos), "\n",
      "These must match. Otherwise the aligned file was made from a different filtered set."
    )
  )
}


aligned_fixed <- data.frame(
  barcode = native1_pos$barcode,
  x = aligned_xy$x,
  y = aligned_xy$y,
  stringsAsFactors = FALSE
)

cat("Preview of corrected file:\n")
print(head(aligned_fixed))


write.csv(
  aligned_fixed,
  output_path,
  row.names = FALSE,
  quote = FALSE
)

cat("Saved corrected file to:\n")
cat(output_path, "\n")