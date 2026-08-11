# Ambient RNA removal with SoupX: automatic contamination estimation.
#
# Builds a SoupChannel from the raw and filtered Cell Ranger matrices, defines
# clusters for soup estimation with a standard Seurat workflow, estimates the
# contamination fraction with autoEstCont, and writes the corrected integer
# counts as a Matrix Market file.
#
# This is the default path, applied to all samples except Ewing sarcoma
# (see soupx_fixed_cf.R).
#
# Usage:
#   Rscript soupx_auto.R <cellranger_outs> <sample_name> <output_dir>
#
#   cellranger_outs  Cell Ranger outs directory for this sample, containing
#                    filtered_feature_bc_matrix and raw_feature_bc_matrix.
#   sample_name      Sample identifier, used as a prefix for the output file.
#   output_dir       Directory the corrected matrix is written to. Created if
#                    it does not exist.
#
# Output:
#   <output_dir>/<sample_name>_soupx.mtx        corrected counts (genes x cells)
#   <output_dir>/contamination_per_sample.txt   estimated contamination fraction
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

### 1. Read the filtered (toc) and raw (tod) count matrices
toc <- Read10X(file.path(sample_path, "filtered_feature_bc_matrix"))
tod <- Read10X(file.path(sample_path, "raw_feature_bc_matrix"))

### 2. Cluster the filtered cells; SoupX uses these groups to estimate the soup
get_soup_groups <- function(sobj) {
  sobj <- CreateSeuratObject(sobj)
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  return(sobj@meta.data[['seurat_clusters']])
}

### 3. Estimate contamination and correct the counts
sc <- SoupChannel(tod, toc)
sc <- setClusters(sc, get_soup_groups(toc))
sc <- autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE)
out <- adjustCounts(sc) 

print(paste0('Sum of counts before: ', sum(toc)))
print(paste0('Sum of counts after: ', sum(out)))
print(paste0('Pct of cells left: ', round(sum(out) / sum(toc), digits = 3)))

### 4. Record the estimated contamination fraction and write the corrected counts
contamination_file <- file.path(out_path, "contamination_per_sample.txt")
write(paste(c(sample_name, sc$fit$rhoEst), collapse = '\t'), file = contamination_file, append = TRUE)

writeMM(t(out), paste0(out_path, '/', sample_name, '_soupx.mtx'))
