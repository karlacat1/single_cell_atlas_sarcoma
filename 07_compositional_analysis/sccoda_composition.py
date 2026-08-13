"""
Compositional analysis with scCODA.

Sum contrast coding expresses each group effect relative to the grand mean across
all groups, but structurally drops one level. The analysis is therefore run twice
with the group ordering rotated so that a different level is dropped each time,
keeping the reference cell type fixed, and the two results are merged to obtain an
estimate for every group.

Input: .h5ad file whose .obs holds the sample identifier, the cell type or cell
state annotation and the grouping covariate.

Output: <prefix>_log2fc.csv, <prefix>_credible.csv and <prefix>_heatmap.pdf
in --output-dir.

Usage:
  python sccoda_composition.py --adata myeloid_final.h5ad --output-dir results/ \
      --group Entity --cell-type-col cell_states \
      --group-order Ewing aRMS eRMS Osteo NRSTS --prefix myeloid
"""

import argparse
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pertpy as pt
import scanpy as sc


def clean_level(name):
    """Strip the sum contrast wrapper, e.g. 'C(Entity, Sum)[S.Ewing]' -> 'Ewing'."""
    m = re.search(r'[ST]\.([^\]]+)', str(name))
    return m.group(1) if m else str(name)


def run_sccoda(adata, group, group_order, cell_type_col, sample_col,
               reference_cell_type, fdr, rng_key=1):
    """Fit one scCODA model with the given group ordering; the last level is dropped."""
    adata = adata.copy()
    adata.obs[group] = pd.Categorical(adata.obs[group], categories=group_order, ordered=True)

    model = pt.tl.Sccoda()
    data = model.load(adata, type='cell_level', generate_sample_level=True,
                      cell_type_identifier=cell_type_col,
                      sample_identifier=sample_col, covariate_obs=[group])
    data = model.prepare(data, modality_key='coda', formula=f'C({group}, Sum)',
                         reference_cell_type=reference_cell_type)
    model.run_nuts(data, modality_key='coda', rng_key=rng_key)
    model.set_fdr(data, est_fdr=fdr)

    effects = model.get_effect_df(data, modality_key='coda')
    credible = model.credible_effects(data, modality_key='coda')

    log2fc = effects['log2-fold change'].unstack(level=0)
    credible = credible.unstack(level=0).reindex(log2fc.index)
    log2fc.columns = [clean_level(c) for c in log2fc.columns]
    credible.columns = [clean_level(c) for c in credible.columns]

    used_reference = data['coda'].uns['scCODA_params']['reference_cell_type']
    return log2fc, credible, used_reference


def plot_heatmap(log2fc, credible, output_path, col_order=None, row_order=None,
                 cmap='RdBu_r', vmax=None, figsize=(5, 5)):
    """log2 fold change heatmap with credible effects marked."""
    if row_order is not None:
        log2fc, credible = log2fc.loc[row_order], credible.loc[row_order]
    if col_order is not None:
        log2fc, credible = log2fc[col_order], credible[col_order]

    if vmax is None:
        vmax = np.nanmax(np.abs(log2fc.values))

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(log2fc.values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect='auto')

    ax.set_xticks(range(len(log2fc.columns)))
    ax.set_xticklabels(log2fc.columns, rotation=50, fontsize=9)
    ax.set_yticks(range(len(log2fc.index)))
    ax.set_yticklabels(log2fc.index, fontsize=9)
    ax.grid(False)

    for i in range(log2fc.shape[0]):
        for j in range(log2fc.shape[1]):
            if bool(credible.values[i, j]):
                ax.text(j, i, '**', ha='center', va='center',
                        fontsize=10, fontweight='bold', color='black')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('log2(observed / expected proportion)')
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Compositional analysis with scCODA.')
    parser.add_argument('--adata', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--prefix', default='sccoda')
    parser.add_argument('--group', required=True,
                        help='Covariate to test, for example Entity or archetype')
    parser.add_argument('--cell-type-col', required=True,
                        help='Column holding the cell type or cell state annotation')
    parser.add_argument('--sample-col', default='sample')
    parser.add_argument('--group-order', nargs='+', required=True,
                        help='Group levels; the last is dropped in the first run')
    parser.add_argument('--reference-cell-type', default='automatic',
                        help='Reference cell type, or automatic to let scCODA choose')
    parser.add_argument('--fdr', type=float, default=0.2)
    parser.add_argument('--col-order', nargs='*', default=None,
                        help='Column order in the heatmap')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    print('-> Number of cells (rows): {} || Number of genes (columns): {}'.format(
        adata.shape[0], adata.shape[1]))

    # run 1: last level of --group-order is dropped
    print('-> Run 1')
    log2fc, credible, reference = run_sccoda(
        adata, args.group, args.group_order, args.cell_type_col, args.sample_col,
        args.reference_cell_type, args.fdr)
    print(f'   reference cell type: {reference}')

    # run 2: rotate the ordering so a different level is dropped, same reference
    rotated = args.group_order[1:] + args.group_order[:1]
    print('-> Run 2')
    log2fc_2, credible_2, _ = run_sccoda(
        adata, args.group, rotated, args.cell_type_col, args.sample_col,
        reference, args.fdr)

    # take the level missing from run 1 out of run 2
    missing = [g for g in args.group_order if g not in log2fc.columns]
    for level in missing:
        log2fc[level] = log2fc_2[level]
        credible[level] = credible_2[level]
    print(f'   recovered from run 2: {", ".join(missing)}')

    log2fc.to_csv(os.path.join(args.output_dir, f'{args.prefix}_log2fc.csv'))
    credible.to_csv(os.path.join(args.output_dir, f'{args.prefix}_credible.csv'))

    plot_heatmap(log2fc, credible.astype(bool),
                 os.path.join(args.output_dir, f'{args.prefix}_heatmap.pdf'),
                 col_order=args.col_order)

    print(f'Done. Results saved in {args.output_dir}')


if __name__ == '__main__':
    main()
