"""
Infer ligand-receptor interactions with LIANA.

Runs LIANA's rank_aggregate method, which combines the predictions of several
individual scoring methods into consensus magnitude and specificity ranks by
robust rank aggregation. For both ranks, lower values indicate more relevant
interactions.

The analysis can be run on the whole cohort at once, or separately within each
level of a grouping variable: per sample, per sarcoma entity, or per TME
archetype. Only --split-by changes between these; everything else is identical.

Ligands or receptors expressed in fewer than 10% of cells per cluster are
excluded (expr_prop, LIANA default).

Input:
  An .h5ad file holding log-normalized counts over the full gene set. The
  notebook this script is derived from passed `adata.raw.to_adata()`; use
  --use-raw to reproduce that behaviour from an object where the full gene set
  is stored in .raw.

Output (written to --output-dir):
  Whole cohort:   uns_liana_<groupby>.pkl
  Split analyses: <level>/uns_liana_<level>_<groupby>.pkl

  Each pickle holds the LIANA result DataFrame, with columns including source,
  target, ligand_complex, receptor_complex, magnitude_rank and specificity_rank.
  These files are the input to interaction_grouping.py.

Usage:
  # per sample
  python run_liana.py --adata atlas.h5ad --output-dir ccc/per_sample/ \
      --groupby cell_types --split-by sample --use-raw

  # per sarcoma entity
  python run_liana.py --adata atlas.h5ad --output-dir ccc/per_entity/ \
      --groupby cell_types --split-by entity --use-raw

  # whole cohort
  python run_liana.py --adata atlas.h5ad --output-dir ccc/all/ \
      --groupby cell_types --use-raw
"""

import argparse
import os
import pickle

import scanpy as sc
from liana.method import rank_aggregate


def run_liana(adata, groupby, key_added, output_path):
    """
    Run rank_aggregate on one dataset and pickle the resulting table.

    :param adata: AnnData holding log-normalized counts
    :param groupby: obs column with the cell type or cell state annotation
    :param key_added: key under which LIANA stores its results in adata.uns
    :param output_path: path of the pickle file to write
    """
    rank_aggregate(adata, groupby=groupby, return_all_lrs=False,
                   use_raw=False, verbose=True, key_added=key_added)

    with open(output_path, 'wb') as f:
        pickle.dump(adata.uns[key_added], f)

    print(f'   {adata.uns[key_added].shape[0]} interactions written to {output_path}')


def main():
    parser = argparse.ArgumentParser(
        description='Infer ligand-receptor interactions with LIANA.')
    parser.add_argument('--adata', required=True,
                        help='Path to the .h5ad file')
    parser.add_argument('--output-dir', required=True,
                        help='Directory for the LIANA result pickles')
    parser.add_argument('--groupby', default='cell_types',
                        help='Column in adata.obs holding the cell type or cell state '
                             'annotation to compute interactions between')
    parser.add_argument('--split-by', default=None,
                        help='Column in adata.obs to split the analysis by, for example '
                             'sample, entity or archetype. Omit to run on the whole cohort.')
    parser.add_argument('--use-raw', action='store_true',
                        help='Take the full gene set from adata.raw before running LIANA')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print('-> Reading dataset')
    adata = sc.read_h5ad(args.adata)
    if args.use_raw:
        adata = adata.raw.to_adata()
    print('   Number of cells (rows): {} || Number of genes (columns): {}'.format(
        adata.shape[0], adata.shape[1]))

    if args.groupby not in adata.obs.columns:
        raise ValueError(f'--groupby column not found in adata.obs: {args.groupby}')

    # whole cohort
    if args.split_by is None:
        key_added = f'liana_{args.groupby}'
        output_path = os.path.join(args.output_dir, f'uns_{key_added}.pkl')
        run_liana(adata, args.groupby, key_added, output_path)
        print('Done.')
        return

    # split analysis
    if args.split_by not in adata.obs.columns:
        raise ValueError(f'--split-by column not found in adata.obs: {args.split_by}')

    levels = adata.obs[args.split_by].unique()
    print(f'-> Running LIANA per {args.split_by} ({len(levels)} levels)')

    failed = []
    for level in levels:
        print(level)
        adata_subgroup = adata[adata.obs[args.split_by] == level].copy()

        level_dir = os.path.join(args.output_dir, str(level))
        os.makedirs(level_dir, exist_ok=True)

        key_added = f'liana_{level}_{args.groupby}'
        output_path = os.path.join(level_dir, f'uns_{key_added}.pkl')

        try:
            run_liana(adata_subgroup, args.groupby, key_added, output_path)
        except Exception as e:
            # most often because a level has too few cells, or too few cell types,
            # for the aggregate methods to run
            print(f'   CCC not possible for {level}: {e}')
            failed.append(str(level))

    if failed:
        print(f'\nNo results for {len(failed)} of {len(levels)} levels: {", ".join(failed)}')
    print('Done.')


if __name__ == '__main__':
    main()
