# Environments

Two Python environments and one R installation were used.

## Python

Both environments run **Python 3.11.6**. They are separate because scCODA and
pertpy pin a different NumPy major version than the rest of the analysis.

| File | Environment | Used by |
|------|-------------|---------|
| `sc_new_requirements.txt` | `sc_new` | everything except the compositional analysis |
| `sccoda_env_requirements.txt` | `sccoda_env` | `07_compositional_analysis/sccoda_composition.py` |

To recreate either:

```bash
conda create -n sc_new python=3.11.6
conda activate sc_new
pip install -r envs/sc_new_requirements.txt
```

Key package versions in `sc_new`:

| Package | Version |
|---------|---------|
| scanpy | 1.11.5 |
| anndata | 0.12.10 |
| scrublet | 0.2.3 |
| scvi-tools | 1.4.2 |
| torch | 2.10.0 |
| harmonypy | 0.0.9 |
| leidenalg | 0.11.0 |
| umap-learn | 0.5.11 |
| celltypist | 1.6.2 |
| gseapy | 1.0.3 |
| liana | 1.7.1 |
| omnipath | 1.0.12 |
| scib | 1.1.7 |
| lifelines | 0.30.3 |
| scipy | 1.17.1 |
| numpy | 1.26.4 |
| pandas | 2.2.2 |
| scikit-learn | 1.8.0 |
| statsmodels | 0.14.6 |

Key package versions in `sccoda_env`:

| Package | Version |
|---------|---------|
| scCODA | 0.1.9 |
| pertpy | 1.0.3 |
| scanpy | 1.11.5 |
| numpy | 2.4.3 |
| pandas | 2.3.3 |

The requirements files list the packages present in each environment, not a
minimal dependency set. Packages installed by conda rather than pip are listed
as comments at the top of each file and are resolved automatically by pip.

## R

**R 4.4.0** (module `R/4.4.0-GCCcore-14.1.0`; the analysis also runs under 4.4.x).

| Package | Version | Used by |
|---------|---------|---------|
| SoupX | 1.6.2 | `02_preprocessing/soupx/` |
| infercnv | 1.20.0 | `03_cnv_analysis/run_infercnv.R` |
| SCEVAN | 1.0.1 | `03_cnv_analysis/run_scevan.R` |
| Seurat | — | `02_preprocessing/soupx/`, `03_cnv_analysis/` |
| Matrix | 1.7-3 | |

```r
install.packages(c("Seurat", "Matrix"))
remotes::install_github("constantAmateur/SoupX")
BiocManager::install("infercnv")
remotes::install_github("AntonioDeFalco/SCEVAN")
```

## Other software

| Software | Version |
|----------|---------|
| Cell Ranger | 6.1.1 |
| MACS iQ View | 1.3.2 |

## Hardware

scVI and scANVI training in `05_integration/run_scanvi.py` used a GPU. The
remaining steps run on CPU. inferCNV and SCEVAN are memory-intensive and were
run on a compute cluster.
