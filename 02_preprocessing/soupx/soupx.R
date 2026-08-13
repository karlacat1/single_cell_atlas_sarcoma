library(SoupX)
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("Usage: Rscript soupx.R <cellranger_outs> <sample_name> <output_dir> --auto\n       Rscript soupx.R <cellranger_outs> <sample_name> <output_dir> --fixed <cf>", call. = FALSE)
}

sample_path <- args[1]
sample_name <- args[2]
out_path <- args[3]
mode <- args[4]

if (!mode %in% c("--auto", "--fixed")) {
    stop("Mode must be --auto or --fixed.", call. = FALSE)
}

if (mode == "--fixed") {
    if (length(args) < 5) stop("--fixed requires a contamination fraction, e.g. --fixed 0.4", call. = FALSE)
    cf <- as.numeric(args[5])
    if (!is.finite(cf) || cf <= 0 || cf >= 1) stop("Contamination fraction must be > 0 and < 1.", call. = FALSE)
}

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

print(paste0("sample path: ", sample_path))
print(paste0("sample name: ", sample_name))
print(paste0("output path: ", out_path))
print(paste0("mode: ", mode))

if (mode == "--fixed") print(paste0("contamination fraction: ", cf))


### 1. Read Cell Ranger matrices

toc <- Read10X(file.path(sample_path, "filtered_feature_bc_matrix"))
tod <- Read10X(file.path(sample_path, "raw_feature_bc_matrix"))

if (is.list(toc)) toc <- toc[["Gene Expression"]]
if (is.list(tod)) tod <- tod[["Gene Expression"]]

if (any(!is.finite(toc@x))) stop("Filtered matrix contains NA/NaN/Inf values.", call. = FALSE)
if (any(!is.finite(tod@x))) stop("Raw matrix contains NA/NaN/Inf values.", call. = FALSE)


### 2. Cluster filtered cells

get_soup_groups <- function(sobj) {
    sobj <- CreateSeuratObject(sobj)
    sobj <- NormalizeData(sobj, verbose = FALSE)
    sobj <- FindVariableFeatures(sobj, nfeatures = min(2000, nrow(sobj)), verbose = FALSE, selection.method = "vst")
    sobj <- ScaleData(sobj, verbose = FALSE)

    npcs <- min(20, nrow(sobj) - 1, ncol(sobj) - 1)

    if (npcs < 2) {
        clusters <- rep("1", ncol(sobj))
        names(clusters) <- colnames(sobj)
        return(clusters)
    }

    sobj <- RunPCA(sobj, npcs = npcs, verbose = FALSE)
    sobj <- FindNeighbors(sobj, dims = 1:npcs, verbose = FALSE)
    sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)

    clusters <- sobj@meta.data[["seurat_clusters"]]
    names(clusters) <- rownames(sobj@meta.data)

    return(clusters)
}

clusters <- get_soup_groups(toc)

if (anyNA(clusters)) stop("Cluster assignments contain NA values.", call. = FALSE)

print(paste0("number of cells: ", ncol(toc)))
print(paste0("number of clusters: ", length(unique(clusters))))


### 3. Create SoupChannel and set clusters

sc <- SoupChannel(tod, toc)
sc <- setClusters(sc, clusters)


### 4. Estimate/set contamination

if (mode == "--auto") {

    soup_est <- sc$soupProfile$est

    if (length(soup_est) == 0 || any(!is.finite(soup_est))) {
        stop("SoupX produced invalid values in soupProfile$est; autoEstCont() cannot estimate contamination.", call. = FALSE)
    }

    sc <- autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE)

    rho <- sc$fit$rhoEst

    if (length(rho) != 1 || !is.finite(rho) || rho <= 0 || rho >= 1) {
        stop("autoEstCont() returned an invalid contamination fraction.", call. = FALSE)
    }

    print(paste0("estimated contamination: ", rho))

    contamination_file <- file.path(out_path, "contamination_per_sample.txt")
    write(paste(c(sample_name, rho), collapse = "\t"), file = contamination_file, append = TRUE)

    out <- adjustCounts(sc, method = "subtraction", roundToInt = TRUE)
    output_file <- file.path(out_path, paste0(sample_name, "_soupx.mtx"))

} else {

    print(paste0("set contamination to: ", cf, " (", cf * 100, "%)"))

    sc <- setContaminationFraction(sc, cf)

    out <- adjustCounts(sc, method = "subtraction", roundToInt = FALSE)
    output_file <- file.path(out_path, paste0(sample_name, "_soupx_", round(cf * 100), "_CF.mtx"))
}


### 5. Validate and write output

if (any(!is.finite(out@x))) stop("Corrected matrix contains NA/NaN/Inf values.", call. = FALSE)

out@x[out@x < 0] <- 0

print(paste0("Sum of counts before: ", sum(toc)))
print(paste0("Sum of counts after: ", sum(out)))
print(paste0("Pct of cells left: ", round(sum(out) / sum(toc), digits = 3)))

writeMM(t(out), output_file)

print(paste0("output: ", output_file))
