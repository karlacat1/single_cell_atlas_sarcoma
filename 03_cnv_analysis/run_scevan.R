#!/usr/bin/env Rscript
#
# Run SCEVAN on a single 10x Genomics sample.
#
# Infers copy-number alterations and classifies malignant versus non-malignant
# cells. SUBCLONES = FALSE reproduces the original analysis. The resulting
# scevan_results.csv feeds classify_malignant_cells.py.
#
# Usage:
#   Rscript run_scevan.R <sample_matrix> <sample_name> <output_dir> [cores]
#
#   sample_matrix  Directory with the 10x matrix files (matrix.mtx,
#                  barcodes.tsv, features.tsv), as written by the preprocessing
#                  pipeline into matrix_files/.
#   sample_name    Sample identifier, passed to SCEVAN and used in output files.
#   output_dir     Directory for scevan_results.csv.
#   cores          Optional number of CPU cores (default: 15).
#
# Requires: SCEVAN 1.0.1, Seurat (R 4.4.x)

# -------------------------------------------------------------------------
# 1. Check required packages
# -------------------------------------------------------------------------

required_packages <- c("SCEVAN", "Seurat")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "The following required R packages are not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script.",
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# 2. Parse command-line arguments
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 4) {
  stop(
    paste(
      "Usage:",
      "Rscript scevan_analysis_publication_ready.R",
      "<sample_matrix> <sample_name> <output_dir> [cores]"
    ),
    call. = FALSE
  )
}

sample_path <- normalizePath(args[[1]], mustWork = TRUE)
sample_name <- args[[2]]
output_path <- normalizePath(args[[3]], mustWork = FALSE)

# Use 15 CPU cores by default, matching the original analysis.
n_cores <- if (length(args) == 4) {
  suppressWarnings(as.integer(args[[4]]))
} else {
  15L
}

if (is.na(n_cores) || n_cores < 1) {
  stop("<cores> must be a positive integer.", call. = FALSE)
}

if (!nzchar(sample_name)) {
  stop("<sample_name> must not be empty.", call. = FALSE)
}

# Create the output directory if it does not already exist.
if (!dir.exists(output_path)) {
  dir.create(
    output_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# -------------------------------------------------------------------------
# 3. Report analysis parameters
# -------------------------------------------------------------------------

message("Starting SCEVAN analysis")
message("Sample matrix : ", sample_path)
message("Sample name   : ", sample_name)
message("Output folder : ", output_path)
message("CPU cores     : ", n_cores)


# -------------------------------------------------------------------------
# 4. Read the 10x Genomics expression matrix
# -------------------------------------------------------------------------

# Read10X() expects a directory containing the 10x Genomics matrix files,
# typically including matrix.mtx, barcodes.tsv, and features.tsv.
raw_data <- Seurat::Read10X(data.dir = sample_path)


# -------------------------------------------------------------------------
# 5. Run SCEVAN
# -------------------------------------------------------------------------

# SCEVAN is used to infer copy-number alterations and classify malignant
# versus non-malignant cells.
#
# SUBCLONES = FALSE reproduces the original analysis and disables subclone
# inference.
scean_results <- SCEVAN::pipelineCNA(
  raw_data,
  sample = sample_name,
  par_cores = n_cores,
  SUBCLONES = FALSE
)


# -------------------------------------------------------------------------
# 6. Save SCEVAN results
# -------------------------------------------------------------------------

results_file <- file.path(
  output_path,
  "scevan_results.csv"
)

write.csv(
  scean_results,
  file = results_file,
  row.names = TRUE
)

message("SCEVAN analysis completed successfully.")
message("Results written to: ", results_file)
