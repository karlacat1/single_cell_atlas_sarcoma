# Reference resources

The analysis scripts read five external reference files from this directory.
They are not redistributed here because of their licensing terms, with the
exception of the cell cycle gene list.

Scripts locate this directory through the `SARCOMA_RESOURCES` environment
variable, falling back to `../resources` relative to the script:

```bash
export SARCOMA_RESOURCES=/path/to/single_cell_atlas_sarcoma/resources
```

## Required files

| File | Used by | Source |
|------|---------|--------|
| `KEGG_RIBOSOME.v2023.1.Hs.txt` | `preprocessing.py`, `merge_samples.py` | MSigDB, gene set `KEGG_RIBOSOME`, human collection v2023.1 |
| `regev_lab_cell_cycle_genes.txt` | `preprocessing.py`, `merge_samples.py` | Tirosh et al. 2016; included in this directory |
| `PanglaoDB_markers_27_Mar_2020.tsv` | `cell_annotation.py` | PanglaoDB, marker file dated 27 March 2020 |
| `msigdb/c4.all.v2023.1.Hs.symbols.gmt` | `cell_annotation.py` | MSigDB v2023.1, collection C4 (computational gene sets) |
| `msigdb/c5.all.v2023.1.Hs.symbols.gmt` | `cell_annotation.py` | MSigDB v2023.1, collection C5 (ontology gene sets) |
| `msigdb/c6.all.v2023.1.Hs.symbols.gmt` | `cell_annotation.py` | MSigDB v2023.1, collection C6 (oncogenic signatures) |
| `msigdb/c8.all.v2023.1.Hs.symbols.gmt` | `cell_annotation.py` | MSigDB v2023.1, collection C8 (cell type signatures) |
| `gene_order_file.txt` | `run_infercnv.R` | GENCODE v27 (hg38), formatted as gene, chromosome, start, stop |

## Where to obtain them

**MSigDB** (`c4`, `c5`, `c6`, `c8` symbol GMTs and the KEGG ribosome gene set) —
https://www.gsea-msigdb.org/gsea/msigdb. Registration is required and the terms
of use prohibit redistribution. Download the human collections for release
v2023.1 and place the four GMT files in `resources/msigdb/`. The ribosome gene
set is obtained by searching for `KEGG_RIBOSOME` and downloading it in the
grp/txt format, which carries two header lines that the scripts skip.

**PanglaoDB** — https://panglaodb.se/markers.html. Download the marker file and
save it as `PanglaoDB_markers_27_Mar_2020.tsv`. Later versions of the file will
give different annotations from those reported here.

**GENCODE gene order file for inferCNV** — derived from the GENCODE v27 human
annotation, https://www.gencodegenes.org/human/release_27.html. inferCNV expects
a tab-separated file with no header and four columns: gene name, chromosome,
start, stop. The inferCNV documentation describes how to generate it from a GTF,
https://github.com/broadinstitute/infercnv/wiki.

**Cell cycle genes** — `regev_lab_cell_cycle_genes.txt` can be downloaded here: https://www.dropbox.com/s/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip?dl=1. 

## Directory layout

```
resources/
├── README.md
├── regev_lab_cell_cycle_genes.txt
├── KEGG_RIBOSOME.v2023.1.Hs.txt
├── PanglaoDB_markers_27_Mar_2020.tsv
├── gene_order_file.txt
└── msigdb/
    ├── c4.all.v2023.1.Hs.symbols.gmt
    ├── c5.all.v2023.1.Hs.symbols.gmt
    ├── c6.all.v2023.1.Hs.symbols.gmt
    └── c8.all.v2023.1.Hs.symbols.gmt
```
