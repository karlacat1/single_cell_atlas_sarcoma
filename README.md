# A single-cell atlas of cancer-educated ecotypes across high-risk pediatric sarcomas

Analysis code for the single-nucleus RNA sequencing atlas of high-risk pediatric
sarcomas: 89 tumor samples from 62 patients across five diagnosis groups,
profiled by 10x Genomics 5' single-nucleus RNA sequencing and MACSima imaging
cyclic staining.

## Contents

Directories follow the order of the analysis.

| Directory | Contents |
|-----------|----------|
| `01_alignment` | Cell Ranger alignment and count matrix generation |
| `02_preprocessing/soupx` | Ambient RNA removal with SoupX (R) |
| `02_preprocessing/qc_annotation_pipeline` | Per-sample quality control, clustering and automated annotation |
| `03_cnv_analysis` | Copy-number inference with SCEVAN and inferCNV (R), and the CNV-based malignant call |
| `04_cell_type_annotation` | Marker-guided curation of the per-sample annotation |
| `05_integration` | Merging of all samples, and integration with Harmony and scVI/scANVI |
| `06_compartments` | Cell state characterization of the stromal, myeloid and T cell compartments |
| `07_compositional_analysis` | Differential cell type abundance with scCODA |
| `08_archetypes` | Definition of tumor microenvironment archetypes |
| `09_interactions` | Ligand-receptor inference with LIANA, and interaction-based sample grouping |
| `10_survival` | Event-free survival analysis and group transitions over time |
| `envs` | Python and R environment specifications |
| `resources` | Where to obtain the external reference files |
| `demo` | A worked example on a single sample |

## 1. System requirements

**Operating system.** Developed and run on Linux (Debian 12 bookworm, x86_64). The Python
code has no platform-specific dependencies; the R packages are available for
Linux, macOS and Windows.

**Software.** Python 3.11.6 in two environments, R 4.4.0, Cell Ranger 6.1.1,
and MACS iQ View 1.3.2 for the spatial proteomic data. Package versions are
pinned in `envs/`, which also lists the key versions in full. There are no
non-standard hardware requirements for most steps; scVI and scANVI training in
`05_integration/run_scanvi.py` used a GPU, and inferCNV and SCEVAN are
memory-intensive and were run on a compute cluster.

**Tested with.** The versions given in `envs/sc_new_requirements.txt` and
`envs/sccoda_env_requirements.txt`.

## 2. Installation guide

```bash
git clone https://github.com/karlacat1/single_cell_atlas_sarcoma.git
cd single_cell_atlas_sarcoma

conda create -n sc_new python=3.11.6
conda activate sc_new
pip install -r envs/sc_new_requirements.txt
```

The compositional analysis uses a second environment, since scCODA and pertpy
pin a different NumPy major version:

```bash
conda create -n sccoda_env python=3.11.6
conda activate sccoda_env
pip install -r envs/sccoda_env_requirements.txt
```

R packages are installed as described in `envs/README.md`.

Typical install time on a normal desktop computer: 10-20 minutes per Python
environment, depending on network speed. Cell Ranger is installed separately
following the 10x Genomics instructions.

Several scripts read external reference files (marker databases, gene sets, a
gene order file). They are not redistributed here; see `resources/README.md` for
what is needed and where to obtain it. Scripts locate them through the
`SARCOMA_RESOURCES` environment variable.

## 3. Demo

`demo/README.md` describes how to run the per-sample preprocessing and
annotation pipeline on a single sample, with the expected output and run time.
The demo takes approximately 10 minutes on one sample of roughly 9,800 nuclei.

The demo data derive from patient material and are subject to controlled access,
and are therefore provided separately rather than in this repository.

## 4. Instructions for use

The analysis runs in the order of the numbered directories. Each script takes
its inputs and outputs as command-line arguments; run any of them with `--help`
for the full list.

**Per sample**

```bash
# 1. Alignment
01_alignment/cellranger_command.sh <sample_id> <fastq_dir> <fastq_prefixes> <output_dir>

# 2. Ambient RNA removal
#    soupx_auto.R for all entities except Ewing sarcoma, where the contamination
#    fraction is set manually with soupx_fixed_cf.R
02_preprocessing/soupx/run_soupx.sh <sample_list> <data_dir> <output_dir> <soupx_script>

# 3. Quality control, clustering and automated annotation
python 02_preprocessing/qc_annotation_pipeline/main.py -d <matrix_dir> -n <sample_id> -o <output_dir>

# 4. Copy-number inference and the CNV-based malignant call
Rscript 03_cnv_analysis/run_scevan.R <matrix_dir> <sample_id> <output_dir>
python 03_cnv_analysis/classify_malignant_cells.py --input-h5ad <h5ad> \
    --scevan-results <csv> --output-h5ad <h5ad> --output-annotations <txt>
Rscript 03_cnv_analysis/run_infercnv.R <sample_dir> <sample_id> <gene_order_file> <annotation_file>

# 5. Marker-guided curation
python 04_cell_type_annotation/annotate_sample.py --adata <h5ad> --config <json> \
    --sample-name <sample_id> --output-dir <dir> --stage explore
python 04_cell_type_annotation/annotate_sample.py --adata <h5ad> --config <json> \
    --sample-name <sample_id> --output-dir <dir> --stage annotate --output-h5ad <h5ad>
```

**Across the cohort**

```bash
# 6. Merge all samples and integrate
python 05_integration/merge_samples.py -d <data_dir> -s <sample_list> -n <name> -m <metadata>
python 05_integration/run_scanvi.py --input-h5ad <h5ad> --output-h5ad <h5ad> ...
python 05_integration/run_harmony.py --input-h5ad <h5ad> --output-h5ad <h5ad> ...

# 7. Cell states within each compartment
python 06_compartments/caf_states.py --adata <h5ad> --output-dir <dir> --output-h5ad <h5ad>
python 06_compartments/macrophage_states.py --adata <h5ad> --output-dir <dir> --output-h5ad <h5ad>
python 06_compartments/tcell_states.py --adata <h5ad> --compartment CD4 --output-dir <dir> --output-h5ad <h5ad>

# 8. Composition, archetypes, interactions and survival
python 07_compositional_analysis/sccoda_composition.py --adata <h5ad> --group Entity ...
python 08_archetypes/archetype_clustering.py --adata <h5ad> --output-dir <dir> --n-clusters 12
python 09_interactions/cell_to_cell_interaction_liana.py --adata <h5ad> --output-dir <dir> --split-by sample
python 09_interactions/interaction_grouping.py --adata <h5ad> --liana-dir <dir> --output-dir <dir>
python 10_survival/survival_km.py --data <csv> --output-dir <dir> --groups C2 C5 C6 C9 C10 C11
```

### Notes on the workflow

The per-sample pipeline in `02_preprocessing/qc_annotation_pipeline` produces an
automated consensus annotation that served as a starting point for manual
review. Malignant status was then assessed against the SCEVAN and inferCNV
profiles together with entity-defining marker expression, and the per-sample
decisions are recorded in the configuration files in `04_cell_type_annotation`.
Final cell type and cell state labels were assigned after cohort-wide
integration, so the labels in the per-sample objects are not those reported in
the manuscript.

Ambient RNA contamination was identified during initial analysis, after which
SoupX correction was applied to all samples and the downstream analysis re-run
on the corrected counts.

Cell state annotation within each compartment was performed by manual review of
the marker scores, enrichment results and expression plots generated by the
scripts in `06_compartments`. Those scripts regenerate the evidence; the
assignment of clusters to cell states was a manual step and is not reproduced in
code.

Cluster numbering depends on random seeds and on the number of principal
components, so a rerun will not necessarily reproduce the original cluster
identities.

## Figures

| Figure | Produced by |
|--------|-------------|
| 1c-g | `05_integration/` |
| 2a-d | `05_integration/`, `06_compartments/` |
| 2e | `07_compositional_analysis/sccoda_composition.py` |
| 3a-g | `06_compartments/caf_states.py` |
| 3h | `07_compositional_analysis/sccoda_composition.py` |
| 4a-h | `06_compartments/tcell_states.py` |
| 4i | `07_compositional_analysis/sccoda_composition.py` |
| 5a-g | `06_compartments/macrophage_states.py` |
| 5h | `07_compositional_analysis/sccoda_composition.py` |
| 6a-b | `08_archetypes/archetype_clustering.py` |
| 6c | `10_survival/survival_km.py` |
| 6d-f | `07_compositional_analysis/sccoda_composition.py` |
| 7a-f | `09_interactions/interaction_grouping.py` |
| 7g | `10_survival/sankey_transitions.py` |
| 7h | `10_survival/survival_km.py` |

Figure 6g shows MACSima imaging cyclic staining data, processed in MACS iQ View;
the processed images and analysis workflows are deposited at Zenodo
(https://doi.org/10.5281/zenodo.21292448).

## Data availability

Sequencing data are available under controlled access; see the Data availability
statement of the manuscript.

## License

MIT. See `LICENSE`.
