"""
Merge the per-sample datasets into a single cohort object.

For each sample, reads the Cell Ranger counts, the SoupX-corrected counts and
the per-sample annotated object, sets the SoupX-corrected matrix as X, applies
uniform quality control (doublet detection, MAD-based outlier filtering), scores
the cell cycle, attaches the clinical metadata, and concatenates all samples.

Ewing sarcoma samples use the SoupX output with the fixed contamination fraction
(*_soupx_40_CF.mtx) where present, and the automatic estimate (*_soupx.mtx)
otherwise.

Quality control is applied identically to every sample here, so that the merged
object is filtered consistently regardless of the per-sample parameters used
during the initial exploratory pass.

External resource files are read from the directory given by the
SARCOMA_RESOURCES environment variable (default: ../resources relative to this
file). See resources/README.md.

Input:
  -d  Directory containing the aligned samples, one <sample>_aligned/outs
      directory per sample.
  -s  Text file listing the samples to merge, one per line.
  -n  Name of the output file.
  -m  Optional semicolon-separated clinical metadata table, indexed by
      'Project #ID'.

Output:
  <data_dir>/01_Merged_samples/<name>/<name>_outer_merged_data.h5ad

Usage:
  python merge_samples.py -d /data/aligned -s samples.txt -n all_89_samples \
      -m metadata.csv
"""

import argparse
import glob
import os
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
from scipy.sparse import csr_matrix
from scipy.stats import median_abs_deviation as mad

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
sc.settings.verbosity = 0

RESOURCES_DIR = os.environ.get(
    'SARCOMA_RESOURCES',
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'resources')
)
RIBO_GENES_FILE = os.path.join(RESOURCES_DIR, 'KEGG_RIBOSOME.v2023.1.Hs.txt')
CELL_CYCLE_GENES_FILE = os.path.join(RESOURCES_DIR, 'regev_lab_cell_cycle_genes.txt')

# obs columns from the exploratory per-sample pass that are not carried forward
DROP_OBS_COLUMNS = [
    'total_counts_mt', 'log1p_total_counts_mt', 'total_counts_ribo',
    'log1p_total_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb',
    'outlier_rib', 'cluster_to_de', 'healthy_vs_tumor_DE', 'over_clustering',
    'conf_score', 'cells', 'highest_confidence', 'cluster_to_tumor',
    'ewing_genes', 'endothelial', 'stromal', 'immune_cells',
    'ewing_literature', 'old_leiden',
]


def sample_path_prep(sample_list, directory_path):
    """
    Locate, for each sample, the raw counts, the SoupX-corrected counts and the
    annotated per-sample object.

    :return: (annotated object paths, raw count paths, SoupX matrix paths)
    """
    sample_dirs = []
    for sample_name in sample_list:
        if not sample_name.endswith(('aligned', 'reanalyzed')):
            sample_path = glob.glob(os.path.join(directory_path, '*', sample_name + '_aligned'))
        else:
            sample_path = glob.glob(os.path.join(directory_path, '*', sample_name))

        if sample_path:
            sample_dirs.append(os.path.join(sample_path[0], 'outs'))
        else:
            print('## No directory found for:', sample_name)

    annotated_paths, raw_counts_paths, soupx_paths = [], [], []

    for sample_dir in sample_dirs:
        raw_counts = os.path.join(sample_dir, 'filtered_feature_bc_matrix')

        # the analysis directory is identified by the cell annotation file it holds
        annotated = glob.glob(os.path.join(sample_dir, '*', '*_cell_annotation.txt'))
        if len(annotated) > 1:
            print('Multiple annotations for:', sample_dir)

        if annotated:
            analysis_dir = os.path.dirname(annotated[0])
        else:
            # fall back to the directory holding the SoupX output
            soupx_any = (glob.glob(os.path.join(sample_dir, '*', '*_soupx_40_CF.mtx'))
                         or glob.glob(os.path.join(sample_dir, '*', '*_soupx.mtx')))
            if not soupx_any:
                print('## No data found for:', sample_dir)
                continue
            analysis_dir = os.path.dirname(soupx_any[0])

        paths = glob.glob(os.path.join(analysis_dir, '*_malignant_vs_non_malignant.h5ad'))
        # fixed contamination fraction where present, automatic estimate otherwise
        soupx_data = (glob.glob(os.path.join(analysis_dir, '*_soupx_40_CF.mtx'))
                      or glob.glob(os.path.join(analysis_dir, '*_soupx.mtx')))

        if paths and soupx_data:
            if len(paths) > 1:
                print('## More than one possibility for sample', analysis_dir)
            annotated_paths.append(paths[0])
            raw_counts_paths.append(raw_counts)
            soupx_paths.append(soupx_data[0])
        else:
            print('## No data found for:', analysis_dir)

    return annotated_paths, raw_counts_paths, soupx_paths


def mad_outlier(adata, metric, nmads, upper_only=False, mt_floor=5):
    """
    Flag cells deviating by more than nmads median absolute deviations from the
    sample median.

    Mitochondrial percentage is treated as upper-tail only. Because many nuclei
    have almost no mitochondrial counts, the MAD-derived threshold can fall very
    low; mt_floor is used as a floor so that cells below this percentage are
    never flagged on this criterion alone.
    """
    M = adata.obs[metric]

    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    threshold = np.median(M) + nmads * mad(M)
    if threshold < mt_floor:
        print('VERY LOW MADS threshold:', threshold)
        threshold = mt_floor

    return M > threshold


def outlier_detection(adata):
    """Initial filtering, QC metric calculation and MAD-based outlier flagging."""
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    ribo_genes = pd.read_table(RIBO_GENES_FILE, skiprows=2, header=None)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['hb'] = adata.var_names.str.contains("^HB[^(P)]")
    adata.var['malat'] = adata.var_names.str.startswith('MALAT')
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                               inplace=True, percent_top=[20], log1p=True)

    adata.obs = adata.obs[[x for x in adata.obs.columns if x not in DROP_OBS_COLUMNS]]

    adata = adata[adata.obs.pct_counts_mt < 8]
    adata.obs['mad_outlier'] = (
        mad_outlier(adata, 'log1p_total_counts', 5)
        | mad_outlier(adata, 'log1p_n_genes_by_counts', 5)
        | mad_outlier(adata, 'pct_counts_in_top_20_genes', 5)
        | mad_outlier(adata, 'pct_counts_mt', 5, upper_only=True)
    )

    return adata


def filter_cells(adata, mads=True):
    """Remove MAD outliers and predicted doublets."""
    print('-> Number of cells (rows): {} || Number of genes (columns): {}'.format(
        adata.shape[0], adata.shape[1]))
    if mads:
        adata = adata[~adata.obs.mad_outlier]
    try:
        adata = adata[~adata.obs['predicted_doublets'], :]
    except (KeyError, TypeError):
        adata.obs['predicted_doublets'] = False
    print('-> Number of cells (rows): {} || Number of genes (columns): {}'.format(
        adata.shape[0], adata.shape[1]))
    return adata


def doublet_detection(adata):
    """Score and flag doublets with Scrublet."""
    adata.var_names_make_unique()
    print('---> Doublet detection:')
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.07)
    out = scrub.scrub_doublets()
    dfg = pd.DataFrame({'doublet_score': out[0], 'predicted_doublets': out[1]},
                       index=adata.obs.index)
    print(dfg.predicted_doublets.sum(), ' predicted doublets')
    adata.obs['doublet_scores'] = dfg['doublet_score']
    adata.obs['predicted_doublets'] = dfg['predicted_doublets']
    return adata


def soupx_contaminating_genes(adata):
    """
    Count, per gene, the number of cells in which SoupX removed counts, by
    comparing the binarized matrices before and after correction.
    """
    before = adata.layers['raw_counts'].copy()
    after = adata.layers['soupX_counts'].copy()

    before.data = np.where(before.data > 0, 1, 0)
    after.data = np.where(after.data > 0, 1, 0)
    changed = ((before - after) == 1).astype(int)
    adata.layers['changed'] = changed
    adata.var['soupX_removed'] = adata.layers['changed'].sum(axis=0).A1

    return adata.var


def cell_cycle_score(adata):
    """Score S and G2/M phase using the gene sets of Tirosh et al. (2016)."""
    print('---> Computing Cell cycle score ')
    cell_cycle_genes = [x.strip() for x in open(CELL_CYCLE_GENES_FILE)]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    scaled_data = sc.pp.scale(adata, copy=True)

    sc.tl.score_genes_cell_cycle(scaled_data, s_genes=s_genes, g2m_genes=g2m_genes)

    adata.obs['S_score'] = scaled_data.obs['S_score']
    adata.obs['G2M_score'] = scaled_data.obs['G2M_score']
    adata.obs['phase'] = scaled_data.obs['phase']
    return adata.obs


def attach_metadata(adata, metadata, sample_name):
    """Attach the clinical annotation for one sample."""
    fields = {
        'Stage': 'Disease_timepoint_histology',
        'Diagnosis': 'Diagnosis (Histology)',
        'Entity': 'Entity',
        'Entity_subgroups': 'Entity_subgroups',
        'Sex': 'Sex',
        'Age': 'Age',
        'CNV_bulk': 'CNV status',
        'biopsy_date': 'Biopsy Date',
        'Location': 'Tumor Location',
        'Qbic_ID': 'QBiC Code',
        '10x_chip': '10X chip run #',
    }
    for obs_col, meta_col in fields.items():
        adata.obs[obs_col] = metadata.loc[sample_name, meta_col]
    return adata


def concat_samples(sample_path_list, raw_matrix_path, soupx_path, sample_name_list,
                   output_directory, output_file, metadata=None):
    """Process each sample and concatenate into a single object."""
    adatas = {}
    os.makedirs(output_directory + output_file, exist_ok=True)
    output_directory = output_directory + output_file

    if metadata is not None:
        metadata = metadata.set_index('Project #ID')

    for sample_path, raw_data, sample_name, soupx_data in zip(
            sample_path_list, raw_matrix_path, sample_name_list, soupx_path):
        print('Adding sample ', sample_path)

        print('Reading directory containing raw data ...')
        adatas[sample_name] = sc.read_10x_mtx(raw_data)
        print('Reading annotated file ...')
        adata_annotated = sc.read_h5ad(sample_path)
        print('Reading soupx decontamination file ...')
        adata_soupx = sc.read_mtx(soupx_data)

        if adata_soupx.shape != adatas[sample_name].shape:
            print('Transposing ..')
            adata_soupx = adata_soupx.transpose()

        # transfer the per-sample annotation onto the raw object
        adatas[sample_name].obs['final_annotation'] = 'unclear'
        for column in adata_annotated.obs.columns:
            adatas[sample_name].obs[column] = 'unclear'
            adatas[sample_name].obs.loc[adata_annotated.obs.index, column] = \
                adata_annotated.obs[column].values
        adatas[sample_name].uns = adata_annotated.uns

        # keep both matrices; the SoupX-corrected counts become X
        adatas[sample_name].layers['raw_counts'] = adatas[sample_name].X
        adatas[sample_name].layers['soupX_counts'] = adata_soupx.X
        adatas[sample_name].X = adatas[sample_name].layers['soupX_counts']

        adatas[sample_name] = doublet_detection(adatas[sample_name])
        adatas[sample_name] = outlier_detection(adatas[sample_name])
        adatas[sample_name].var = soupx_contaminating_genes(adatas[sample_name])
        adatas[sample_name] = filter_cells(adatas[sample_name], mads=True)

        adatas[sample_name].obs = cell_cycle_score(adatas[sample_name])
        adatas[sample_name].obs = adatas[sample_name].obs.astype('str')

        if metadata is not None:
            adatas[sample_name] = attach_metadata(adatas[sample_name], metadata, sample_name)

    print('Concatenating all data...')
    adata = ad.concat(adatas, label='sample', index_unique='_', join='outer', fill_value=0)
    print('-> Number of cells (rows): {} || Number of genes (columns): {}'.format(
        adata.shape[0], adata.shape[1]))

    print('Saving concatenated Data...')
    adata.X = csr_matrix(adata.X)
    adata.write(output_directory + output_file + '_outer_merged_data.h5ad')
    print('Done')


def main():
    parser = argparse.ArgumentParser(prog='Sample_merger')
    parser.add_argument('-d', nargs=1, help='Path to directory containing data', required=True)
    parser.add_argument('-s', nargs=1, help='Path to list of samples to be merged', required=True)
    parser.add_argument('-n', nargs=1, help='Name of output file', required=True)
    parser.add_argument('-m', nargs=1, help='Metadata file')

    args = parser.parse_args()
    directory_path = args.d[0]
    sample_list_path = args.s[0]
    output_file_name = args.n[0]
    metadata = None

    if args.m is not None:
        metadata = pd.read_csv(args.m[0], sep=';', encoding='unicode_escape')

    if not os.path.exists(sample_list_path) or not os.path.exists(directory_path):
        parser.error('Path incorrect!')

    with open(sample_list_path, 'r') as file:
        sample_list = file.read().splitlines()

    sample_paths, raw_matrix_path, soupx_path = sample_path_prep(sample_list, directory_path)

    output_directory = os.path.join(directory_path, '01_Merged_samples')
    os.makedirs(output_directory, exist_ok=True)

    concat_samples(sample_paths, raw_matrix_path, soupx_path, sample_list,
                   output_directory, '/' + output_file_name, metadata)

    print('Done. Results saved in directory ' + output_directory + '\n\n')


if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print('Time: {} min'.format((end - start) / 60))
