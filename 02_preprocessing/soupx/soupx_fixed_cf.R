##### SoupX manual correction script ##### written by Anke King
#
# Ambient RNA removal with SoupX using a manually set contamination fraction.
#
# Applied to Ewing sarcoma samples, in which autoEstCont underestimated the
# contamination: Ewing-associated transcripts remained detectable in
# non-malignant cells after automatic correction. A range of contamination
# fractions was tested and evaluated with biological conservation and
# integration metrics (scib); 0.4 removed Ewing-associated genes from the
# non-malignant compartment while retaining their expression in malignant cells.
#
# This script differs from soupx_auto.R in three respects: the contamination
# fraction is set manually rather than estimated; the SoupChannel is built with
# load10X, which uses the Cell Ranger clustering to define soup groups rather
# than a Seurat clustering; and adjustCounts is called without roundToInt, so
# the corrected counts are not rounded to integers.
#
# Usage:
#   Rscript soupx_fixed_cf.R <cellranger_dir> <sample_name> <output_dir>
#
#   cellranger_dir  Cell Ranger outs directory for this sample, readable by
#                   load10X (filtered and raw matrices plus analysis/clustering).
#   sample_name     Sample identifier, used as a prefix for the output file.
#   output_dir      Directory for the corrected counts. Created if absent.
#
# Output:
#   <output_dir>/<sample_name>_soupx_40_CF.mtx   corrected counts (genes x cells)
#
# Requires: SoupX 1.6.2, Seurat, Matrix (R 4.4.x)

library(SoupX)
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Please supply: \n 1. Cell Ranger outs directory\n
       2. Sample name\n
       3. Output directory", call. = FALSE)
}

sample_path <- args[1]
sample_name <- args[2]
out_path <- args[3]

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

print(paste0('sample path: ', sample_path))
print(paste0('sample name: ', sample_name))
print(paste0('output path: ', out_path))

# Load the Cell Ranger output as a SoupChannel
sc <- load10X(sample_path)

# Set the contamination fraction
cf <- 0.4
cf_per <- as.character(cf * 100)
print(paste0("set contamination to: ", cf_per, "%"))

sc <- setContaminationFraction(sc, cf)

# Adjust counts
out <- adjustCounts(sc)

# Write output matrix
writeMM(t(out), paste0(out_path, '/', sample_name, '_soupx_40_CF.mtx'))
