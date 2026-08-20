"""
Merge per-sample single-cell datasets into a single cohort object.

For each sample, reads the raw counts and (optionally) a SoupX-corrected
counts matrix and a pre-annotated object, sets the SoupX-corrected matrix as
X when available, applies uniform quality control (doublet detection,
MAD-based outlier filtering), scores the cell cycle, optionally attaches
clinical/metadata columns, and concatenates all samples into one AnnData.

Sample inputs are supplied explicitly via a sample sheet, rather than being
inferred from directory-naming conventions -- this keeps the script portable
across projects with different folder layouts.

Sample sheet (CSV or TSV, detected from the file extension unless
--sample-sheet-sep is given) with one row per sample and columns:

  sample          Sample identifier (required)
  raw_counts      Path to raw counts: a 10x mtx directory, a 10x .h5 file,
                   or an .h5ad file (required)
  soupx_counts    Path to a SoupX-corrected matrix, .mtx format (optional)
  annotated       Path to a pre-annotated .h5ad whose obs columns should be
                   transferred onto the raw object by matching cell
                   barcodes (optional)

Example sample sheet (samples.tsv):

  sample     raw_counts                              soupx_counts               annotated
  S001       /data/S001/filtered_feature_bc_matrix   /data/S001/soupx.mtx       /data/S001/S001_annotated.h5ad
  S002       /data/S002/filtered_feature_bc_matrix                              /data/S002/S002_annotated.h5ad

Metadata (optional, --metadata) is a separate table indexed by sample name;
every column present is attached to adata.obs for the matching sample,
broadcast to all its cells. No fixed schema is assumed -- whatever columns
are in the file get attached.

Output:
  <output>.h5ad

Usage:
  python merge_samples.py \
      --sample-sheet samples.tsv \
      --ribo-genes resources/KEGG_RIBOSOME.v2023.1.Hs.txt \
      --cell-cycle-genes resources/regev_lab_cell_cycle_genes.txt \
      --output cohort_merged.h5ad \
      [--metadata metadata.csv --metadata-sep , --metadata-index-col sample_id] \
      [--drop-obs-columns col1,col2] \
      [--min-genes 200] [--min-cells 3] [--mt-pct-max 8] [--mt-floor 5] \
      [--nmads 5] [--doublet-rate 0.07]
"""

import argparse
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


def read_counts_matrix(path):
    """Read a counts matrix from a 10x mtx directory, a 10x .h5 file, or an .h5ad file."""
    if os.path.isdir(path):
        return sc.read_10x_mtx(path)
    if path.endswith('.h5ad'):
        return sc.read_h5ad(path)
    if path.endswith('.h5'):
        return sc.read_10x_h5(path)
    raise ValueError(
        f'Could not determine how to read counts matrix at: {path} '
        '(expected a 10x mtx directory, a .h5 file, or an .h5ad file)'
    )


def read_sample_sheet(path, sep=None):
    """Read the sample sheet, auto-detecting the delimiter from the file extension if not given."""
    if sep is None:
        sep = '\t' if path.endswith(('.tsv', '.txt')) else ','
    df = pd.read_csv(path, sep=sep, dtype=str)
    if 'sample' not in df.columns or 'raw_counts' not in df.columns:
        raise ValueError("Sample sheet must contain at least 'sample' and 'raw_counts' columns.")
    for optional_col in ('soupx_counts', 'annotated'):
        if optional_col not in df.columns:
            df[optional_col] = np.nan
    return df


def mad_outlier(adata, metric, nmads, upper_only=False, mt_floor=None):
    """
    Flag cells deviating by more than nmads median absolute deviations from
    the sample median.

    If upper_only is True, only the upper tail is flagged (appropriate for
    metrics like mitochondrial percentage, where many cells sit near zero and
    a two-sided MAD would flag low-percentage cells for being "too low").
    mt_floor, if given, is a minimum threshold below which flagging never
    happens on this criterion alone -- useful when the MAD-derived threshold
    would otherwise fall implausibly low.
    """
    M = adata.obs[metric]

    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    threshold = np.median(M) + nmads * mad(M)
    if mt_floor is not None and threshold < mt_floor:
        print(f'  MAD threshold for {metric} ({threshold:.3f}) below floor, using floor {mt_floor}')
        threshold = mt_floor

    return M > threshold


def outlier_detection(adata, ribo_genes_file, min_genes, min_cells,
                      mt_pct_max, mt_floor, nmads, mt_prefix='MT-',
                      hb_regex='^HB[^(P)]', malat_prefix='MALAT',
                      drop_obs_columns=None):
    """Initial filtering, QC metric calculation and MAD-based outlier flagging."""
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    ribo_genes = pd.read_table(ribo_genes_file, skiprows=2, header=None)
    adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
    adata.var['hb'] = adata.var_names.str.contains(hb_regex)
    adata.var['malat'] = adata.var_names.str.startswith(malat_prefix)
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                               inplace=True, percent_top=[20], log1p=True)

    if drop_obs_columns:
        adata.obs = adata.obs[[c for c in adata.obs.columns if c not in drop_obs_columns]]

    adata = adata[adata.obs.pct_counts_mt < mt_pct_max]
    adata.obs['mad_outlier'] = (
        mad_outlier(adata, 'log1p_total_counts', nmads)
        | mad_outlier(adata, 'log1p_n_genes_by_counts', nmads)
        | mad_outlier(adata, 'pct_counts_in_top_20_genes', nmads)
        | mad_outlier(adata, 'pct_counts_mt', nmads, upper_only=True, mt_floor=mt_floor)
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


def doublet_detection(adata, expected_doublet_rate):
    """Score and flag doublets with Scrublet."""
    adata.var_names_make_unique()
    print('---> Doublet detection:')
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
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


def cell_cycle_score(adata, cell_cycle_genes_file, n_s_genes=43):
    """
    Score S and G2/M phase.

    cell_cycle_genes_file should list S-phase genes followed by G2/M genes,
    one per line, with the first n_s_genes lines being the S-phase set
    (default 43, matching the Tirosh et al. 2016 gene list).
    """
    print('---> Computing cell cycle score')
    cell_cycle_genes = [x.strip() for x in open(cell_cycle_genes_file)]
    s_genes = cell_cycle_genes[:n_s_genes]
    g2m_genes = cell_cycle_genes[n_s_genes:]
    scaled_data = sc.pp.scale(adata, copy=True)

    sc.tl.score_genes_cell_cycle(scaled_data, s_genes=s_genes, g2m_genes=g2m_genes)

    adata.obs['S_score'] = scaled_data.obs['S_score']
    adata.obs['G2M_score'] = scaled_data.obs['G2M_score']
    adata.obs['phase'] = scaled_data.obs['phase']
    return adata.obs


def attach_metadata(adata, metadata, sample_name):
    """Attach every column of the metadata table to adata.obs for this sample."""
    if sample_name not in metadata.index:
        print(f'  Warning: no metadata row found for sample "{sample_name}", skipping.')
        return adata
    row = metadata.loc[sample_name]
    for col in metadata.columns:
        adata.obs[col] = row[col]
    return adata


def process_sample(sample_row, args, metadata=None):
    """Load, QC, and annotate a single sample."""
    sample_name = sample_row['sample']
    print('Processing sample:', sample_name)

    print('  Reading raw counts...')
    adata = read_counts_matrix(sample_row['raw_counts'])

    has_annotation = isinstance(sample_row.get('annotated'), str) and sample_row['annotated']
    if has_annotation:
        print('  Reading annotated file...')
        adata_annotated = sc.read_h5ad(sample_row['annotated'])
        for column in adata_annotated.obs.columns:
            adata.obs[column] = 'unclear'
            common = adata.obs.index.intersection(adata_annotated.obs.index)
            adata.obs.loc[common, column] = adata_annotated.obs.loc[common, column].values
        adata.uns = adata_annotated.uns

    has_soupx = isinstance(sample_row.get('soupx_counts'), str) and sample_row['soupx_counts']
    adata.layers['raw_counts'] = adata.X
    if has_soupx:
        print('  Reading SoupX-corrected counts...')
        adata_soupx = sc.read_mtx(sample_row['soupx_counts'])
        if adata_soupx.shape != adata.shape:
            print('  Transposing SoupX matrix to match raw counts shape...')
            adata_soupx = adata_soupx.transpose()
        adata.layers['soupX_counts'] = adata_soupx.X
        adata.X = adata.layers['soupX_counts']
    # if no SoupX matrix was supplied, X simply stays as the raw counts

    adata = doublet_detection(adata, args.doublet_rate)
    adata = outlier_detection(
        adata,
        ribo_genes_file=args.ribo_genes,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        mt_pct_max=args.mt_pct_max,
        mt_floor=args.mt_floor,
        nmads=args.nmads,
        drop_obs_columns=args.drop_obs_columns,
    )
    if has_soupx:
        adata.var = soupx_contaminating_genes(adata)
    adata = filter_cells(adata, mads=True)

    adata.obs = cell_cycle_score(adata, args.cell_cycle_genes)
    adata.obs = adata.obs.astype('str')

    if metadata is not None:
        adata = attach_metadata(adata, metadata, sample_name)

    return adata


def merge_samples(sample_sheet, args, metadata=None):
    """Process every sample in the sheet and concatenate into one AnnData."""
    adatas = {}
    for _, row in sample_sheet.iterrows():
        adatas[row['sample']] = process_sample(row, args, metadata=metadata)

    print('Concatenating all samples...')
    adata = ad.concat(adatas, label='sample', index_unique='_', join='outer', fill_value=0)
    print('-> Number of cells (rows): {} || Number of genes (columns): {}'.format(
        adata.shape[0], adata.shape[1]))

    adata.X = csr_matrix(adata.X)
    return adata


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merge per-sample single-cell datasets into a single cohort object.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--sample-sheet', required=True,
                        help='CSV/TSV with columns: sample, raw_counts, [soupx_counts], [annotated]')
    parser.add_argument('--sample-sheet-sep', default=None,
                        help='Delimiter for the sample sheet (default: auto-detect from extension)')
    parser.add_argument('--output', required=True, help='Output .h5ad path')
    parser.add_argument('--ribo-genes', required=True,
                        help='Path to a ribosomal gene list file (one gene per line after 2 header rows)')
    parser.add_argument('--cell-cycle-genes', required=True,
                        help='Path to a cell cycle gene list file (S-phase genes then G2/M genes, one per line)')
    parser.add_argument('--metadata', default=None,
                        help='Optional metadata table to attach to adata.obs')
    parser.add_argument('--metadata-sep', default=',',
                        help='Delimiter for the metadata table (default: ",")')
    parser.add_argument('--metadata-index-col', default=None,
                        help='Column in the metadata table matching the "sample" values in the sample sheet '
                             '(default: first column)')
    parser.add_argument('--drop-obs-columns', default=None,
                        help='Comma-separated obs column names to drop before QC (default: none)')
    parser.add_argument('--min-genes', type=int, default=200)
    parser.add_argument('--min-cells', type=int, default=3)
    parser.add_argument('--mt-pct-max', type=float, default=8.0,
                        help='Cells with pct_counts_mt above this are removed outright before MAD flagging')
    parser.add_argument('--mt-floor', type=float, default=5.0,
                        help='Minimum pct_counts_mt MAD threshold (see mad_outlier docstring)')
    parser.add_argument('--nmads', type=float, default=5.0,
                        help='Number of median absolute deviations for outlier flagging')
    parser.add_argument('--doublet-rate', type=float, default=0.07,
                        help='Expected doublet rate passed to Scrublet')

    args = parser.parse_args()
    args.drop_obs_columns = (
        [c.strip() for c in args.drop_obs_columns.split(',')] if args.drop_obs_columns else []
    )
    return args


def main():
    args = parse_args()

    if not os.path.exists(args.sample_sheet):
        raise FileNotFoundError(f'Sample sheet not found: {args.sample_sheet}')

    sample_sheet = read_sample_sheet(args.sample_sheet, sep=args.sample_sheet_sep)

    metadata = None
    if args.metadata is not None:
        metadata = pd.read_csv(args.metadata, sep=args.metadata_sep, encoding='unicode_escape')
        index_col = args.metadata_index_col or metadata.columns[0]
        metadata = metadata.set_index(index_col)

    adata = merge_samples(sample_sheet, args, metadata=metadata)

    print('Saving merged data...')
    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)
    adata.write(args.output)
    print('Done. Results saved to:', args.output)


if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print('Time: {:.2f} min'.format((end - start) / 60))
