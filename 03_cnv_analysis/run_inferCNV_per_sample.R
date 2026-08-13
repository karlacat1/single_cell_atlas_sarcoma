#!/usr/bin/env Rscript
#
# run_infercnv.R
#
# Purpose:
#   Runs inferCNV (https://github.com/broadinstitute/infercnv) on a single-cell
#   RNA-seq sample to estimate copy number variation (CNV) profiles from gene
#   expression data. Reference ("normal") cells are automatically identified
#   from the cell-type annotation file, preferring T cells when available and
#   otherwise falling back to any non-malignant, non-stromal cell type.
#
# Usage:
#   Rscript run_infercnv.R <sample_path> <sample_name> [gene_order_file]
#
# Arguments:
#   sample_path      Path to the directory containing the sample's analyzed
#                     data. This directory must contain:
#                       - matrix_files/              (10x-format count matrix)
#                       - metadata.csv                (cell metadata; must have
#                                                       a column "X" with cell
#                                                       barcodes as rownames,
#                                                       used for a "leiden"
#                                                       cluster column)
#                       - <sample_name>_cell_annotation.txt
#                                                      (tab-separated, no header,
#                                                       rownames = cell barcodes,
#                                                       column V2 = cell type)
#   sample_name      Sample identifier (used to locate the cell annotation
#                     file and to name outputs).
#   gene_order_file  (Optional) Path to the inferCNV gene ordering file
#                     (chromosome/start/stop per gene). If not supplied, the
#                     script looks for a file named "gene_order_file.txt" in
#                     the same directory as this script.
#
# Output:
#   inferCNV results are written to <sample_path>/infercnv_output/
#
# Requirements:
#   R packages: infercnv, Seurat, devtools
#
# Notes:
#   - If a local R library path is needed (e.g. on an HPC cluster without
#     write access to the default library), set the R_LIBS_USER environment
#     variable before running this script, or uncomment and edit the
#     .libPaths() line below.

# .libPaths("/path/to/your/R_packages")  # uncomment and edit if needed

suppressPackageStartupMessages({
  library(devtools)
  library(infercnv)
  library(Seurat)
})

## ---- Parse command-line arguments -----------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript run_infercnv.R <sample_path> <sample_name> [gene_order_file]\n",
    "  sample_path      Path to the directory with the sample's analyzed data\n",
    "  sample_name      Sample identifier (without any '_aligned' suffix)\n",
    "  gene_order_file  Optional path to the inferCNV gene order file\n",
    call. = FALSE
  )
}

sample_path <- args[1]
sample_name <- args[2]

# Default gene order file: look next to this script unless one was supplied.
if (length(args) >= 3) {
  path_to_gene_order_file <- args[3]
} else {
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
  )
  path_to_gene_order_file <- file.path(script_dir, "gene_order_file.txt")
}

cat(sprintf("Sample path: %s\n", sample_path))
cat(sprintf("Sample name: %s\n", sample_name))
cat(sprintf("Gene order file: %s\n", path_to_gene_order_file))

if (!dir.exists(sample_path)) {
  stop(sprintf("Sample path does not exist: %s", sample_path), call. = FALSE)
}
if (!file.exists(path_to_gene_order_file)) {
  stop(sprintf("Gene order file not found: %s", path_to_gene_order_file), call. = FALSE)
}

setwd(sample_path)

## ---- Load input data -------------------------------------------------------

matrix_dir     <- file.path(sample_path, "matrix_files")
metadata_path  <- file.path(sample_path, "metadata.csv")
annotation_path <- file.path(sample_path, paste0(sample_name, "_cell_annotation.txt"))

if (!dir.exists(matrix_dir))      stop(sprintf("Matrix directory not found: %s", matrix_dir), call. = FALSE)
if (!file.exists(metadata_path))  stop(sprintf("Metadata file not found: %s", metadata_path), call. = FALSE)
if (!file.exists(annotation_path)) stop(sprintf("Cell annotation file not found: %s", annotation_path), call. = FALSE)

# Cell metadata (must contain a "X" column with cell barcodes and a "leiden"
# clustering column used below to disambiguate cell subtypes).
metadata_file <- read.csv(metadata_path)
rownames(metadata_file) <- metadata_file$X

# Raw count matrix (10x Genomics format) -> Seurat object.
raw_counts <- Read10X(data.dir = matrix_dir)
sobj <- CreateSeuratObject(counts = raw_counts, meta.data = metadata_file)

# Extract raw counts matrix from the Seurat object for inferCNV.
counts_matrix <- GetAssayData(sobj)

# Cell-type annotation file (tab-separated: cell barcode, cell type).
celltype_annotation_file <- read.csv(
  annotation_path,
  sep = "\t",
  header = FALSE,
  row.names = 1
)

## ---- Determine reference (non-malignant) cell groups -----------------------

# Prefer T cells as the reference group if present in the annotation.
ref_group <- unique(celltype_annotation_file$V2)
ref_group <- ref_group[ref_group == "T cells"]

# Otherwise, fall back to any cell type that is not malignant or a known
# stromal/structural cell type.
non_reference_types <- c(
  "malignant", "Fibroblasts", "Pericytes", "Stromal cells", "Muscle cells",
  "Epithelial cells", "Endothelial cells", "Smooth muscle cells",
  "unknown", "unclear"
)

if (length(ref_group) == 0) {
  ref_group <- unique(celltype_annotation_file$V2)
  ref_group <- ref_group[!(ref_group %in% non_reference_types)]
}

if (length(ref_group) == 0) {
  stop("No suitable reference cell group could be identified from the annotation file.", call. = FALSE)
}

## ---- Disambiguate non-reference cell types by cluster -----------------------

# For any cell type that is not part of the reference group, append its
# Leiden cluster number to the cell type label (e.g. "malignant_3"). This
# lets inferCNV treat distinct clusters of the same broad cell type
# (e.g. different malignant subclones) separately.
metadata_file$celltype <- celltype_annotation_file$V2

for (cell_type in unique(celltype_annotation_file$V2)) {
  if (!(cell_type %in% ref_group)) {
    is_type <- metadata_file$celltype == cell_type
    metadata_file$celltype[is_type] <- paste(
      metadata_file$celltype[is_type],
      metadata_file$leiden[is_type],
      sep = "_"
    )
  }
}

celltype_annotation_file$V2 <- metadata_file$celltype

## ---- Run inferCNV -----------------------------------------------------------

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file   = celltype_annotation_file,
  delim              = "\t",
  gene_order_file    = path_to_gene_order_file,
  ref_group_names    = ref_group
)

output_dir <- file.path(sample_path, "infercnv_output")

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff              = 0.1,   # 1 for Smart-seq2 data; 0.1 for 10x Genomics data
  out_dir             = output_dir,
  min_cells_per_gene  = 3,
  cluster_by_groups   = TRUE,
  denoise             = TRUE,
  num_threads         = 40,
  HMM                 = TRUE
)

cat("Done.\n")
