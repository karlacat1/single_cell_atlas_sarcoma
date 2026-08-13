# Demo

Runs the per-sample preprocessing and annotation pipeline
(`02_preprocessing/qc_annotation_pipeline/`) on a single sarcoma sample.

The demo covers quality control and filtering, doublet removal, normalization,
feature selection, PCA, clustering, UMAP, and the automated consensus cell type
annotation. It does not cover the cohort-level steps (integration, archetypes,
cell-cell interactions, survival), which require the full 89-sample dataset.

## Demo data

The demo input is the Cell Ranger filtered feature-barcode matrix for one
sample, provided separately with the submission. It is not included in this
repository because it derives from patient material and is subject to
controlled access.

Place it as:

```
demo/
├── data/
│   └── filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
└── resources/
    ├── KEGG_RIBOSOME.v2023.1.Hs.txt
    ├── regev_lab_cell_cycle_genes.txt
    ├── PanglaoDB_markers_27_Mar_2020.tsv
    └── msigdb/
```

The reference files in `demo/resources/` are the ones described in
`resources/README.md`; they are not redistributed here and must be downloaded
from their respective sources.

## Prerequisites

1. The `sc_new` environment, created as described in `envs/README.md`.
2. The reference files listed in `resources/README.md`, placed in
   `demo/resources/` and pointed to by `SARCOMA_RESOURCES`.

## Running the demo

```bash
conda activate sc_new
export SARCOMA_RESOURCES=/path/to/single_cell_atlas_sarcoma/demo/resources

cd 02_preprocessing/qc_annotation_pipeline

python main.py \
    -d ../../demo/data/filtered_feature_bc_matrix \
    -n DEMO_SAMPLE \
    -o ../../demo/output
```

To run only the quality control and preprocessing steps, which is faster and
does not require the PanglaoDB and MSigDB files:

```bash
python main.py \
    -d ../../demo/data/filtered_feature_bc_matrix \
    -n DEMO_SAMPLE \
    -o ../../demo/output \
    -pr
```

## Expected output

```
demo/output/
├── DEMO_SAMPLE_non_log.h5ad                     counts before log transformation
├── DEMO_SAMPLE_preprocessed.h5ad                preprocessed, clustered object
├── DEMO_SAMPLE_malignant_vs_non_malignant.h5ad  annotated object
├── DEMO_SAMPLE_malignant_annotation.csv         cluster annotation for inferCNV
├── metadata.csv                                 per-cell metadata
├── matrix_files/                                MEX export for the R scripts
├── 1.preprocessing/
│   ├── DEMO_SAMPLE_top_20_highest_expr.png
│   ├── DEMO_SAMPLE_violin_qc_plots_initial.png
│   ├── DEMO_SAMPLE_violin_qc_plots_after_doublet.png
│   ├── DEMO_SAMPLE_violin_qc_plots_after_filtering.png
│   ├── DEMO_SAMPLE_mt_vs_count_depth.png
│   ├── DEMO_SAMPLE_genes_vs_count_depth.png
│   ├── DEMO_SAMPLE_doublet_compare.png
│   ├── DEMO_SAMPLE_doublet_histogram.png
│   └── DEMO_SAMPLE_pca_elbow.png
└── 2.downstream_analysis/
    ├── DEMO_SAMPLE_leiden_umap2.png
    ├── DEMO_SAMPLE_pca_clustering.png
    ├── DEMO_SAMPLE_cell_cycle_score.png
    ├── DEMO_SAMPLE_DE_genes_wilcoxon.png
    ├── DEMO_SAMPLE_DE_top_marker.csv
    ├── DEMO_SAMPLE_DE_celltype_heatmap.png
    ├── DEMO_SAMPLE_confidence_score.png
    ├── DEMO_SAMPLE_consensus_result.csv
    ├── DEMO_SAMPLE_consensus_result_complete.csv
    ├── DEMO_SAMPLE_count_result.csv
    ├── DEMO_SAMPLE_tumor_marker_expression_umap2.png
    ├── DEMO_SAMPLE_multiple_umap.png
    ├── DEMO_SAMPLE_consensus malignant_umap2.png
    ├── DEMO_SAMPLE_consensus all_umap2.png
    └── gsea/
```

The console reports progress through each step. For the demo sample the run
should report approximately:

| | |
|---|---|
| Cells read | 9,778 |
| Genes read | 36,601 |
| Predicted doublets | 1 |
| Cells after filtering | 9,505 |
| Genes after filtering | 7,202 |
| Leiden clusters | 13 |
| Malignant cells | 3,158 |
| Non-malignant cells | 4,226 |
| Unclear cells | 2,121 |

The same summary is written to `DEMO_SAMPLE_count_result.csv`. Exact numbers
depend on the versions of Scanpy and its dependencies.

## Expected run time

On a cluster worker node, for this sample of roughly 9,800 nuclei:

| Step | Time |
|------|------|
| Preprocessing only (`-pr`) | approximately 2 minutes |
| Full run including annotation | approximately 10 minutes |

The gene set enrichment step accounts for most of the full run: it performs
1,000 permutations against four MSigDB collections for each of the 13 clusters,
at roughly 30 seconds per cluster.

## Interpreting the result

The pipeline produces an automated consensus annotation, combining
marker-gene overlap, CellTypist, PanglaoDB and gene set enrichment. This is a
starting point for manual curation, not the final annotation reported in the
manuscript: malignant status was subsequently reviewed against the inferCNV and
SCEVAN profiles (`03_cnv_analysis/`) and entity-defining marker expression
(`04_cell_type_annotation/`), and the final cell type and cell state labels were
assigned after cohort-wide integration (`05_integration/`, `06_compartments/`).
Cluster numbering depends on the Leiden random seed and on the number of
principal components, so cluster identities in a rerun will not necessarily
match those in the original analysis.

## Troubleshooting

**Segmentation fault, or an OpenBLAS error about too many threads.** On machines
with many cores, OpenBLAS can exceed the number of thread regions it was built
to handle, and the run dies during the differential expression step. Limit the
thread count before starting Python:

```bash
export OPENBLAS_NUM_THREADS=8
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
```

The variables must be set before Python starts, since the numerical libraries
read them at import. If the gene set enrichment step then runs out of threads,
lower `GSEA_THREADS` in `cell_annotation.py` from its default of 25 to match.

**Warnings during the run.** Several messages are expected and do not indicate
a problem: Scanpy reporting `... storing '<column>' as categorical` when the
annotated object is written, CellTypist noting that it will use `.raw.X`
instead of `.X`, gseapy warning about duplicated values in the preranked
statistics, and three cell cycle genes (`MLF1IP`, `FAM64A`, `HN1`) being absent
from the dataset because the reference list uses older gene symbols.

**Very few doublets detected.** Scrublet typically calls only a handful of
doublets in single-nucleus data, and reports an estimated overall rate well
above the expected 7%. Nuclei have lower transcriptional complexity than whole
cells, so simulated doublets separate less cleanly and the automatic threshold
ends up high. This is expected and does not indicate a failed run.

**`ModuleNotFoundError`.** All six pipeline files must be in the same directory,
and the working directory must be that directory, so that Python finds the
modules `main.py` imports.
