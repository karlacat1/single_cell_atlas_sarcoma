#!/usr/bin/env python3
"""
CNV-based malignant / non-malignant call from SCEVAN results.

Each Leiden cluster is assigned the most frequent SCEVAN class. Clusters called
"tumor" or "filtered" are labelled malignant; SCEVAN frequently assigns tumor
cells to "filtered" rather than "tumor", so both are treated as malignant here.

This is a permissive first pass. The final call is made by inspection of the
inferCNV profiles together with entity-defining marker expression, recorded per
sample in annotate_sample.py.

Input: .h5ad with adata.obs["leiden"], and the SCEVAN results CSV indexed by
cell barcode with a "class" column.

Output: .h5ad with adata.obs["scevan_class"] and adata.obs["cnv_annotation"],
plus a tab-separated annotation table. Feeds annotate_sample.py.

Usage:
  python classify_malignant_cells.py --input-h5ad sample.h5ad \
      --scevan-results scevan_results.csv \
      --output-h5ad sample_cnv_classified.h5ad \
      --output-annotations sample_cnv_annotation.txt
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

MALIGNANT_CLASSES = ('tumor', 'filtered')


def add_scevan_predictions(adata, scevan):
    """Add the per-cell SCEVAN class to adata.obs."""
    common = adata.obs.index.intersection(scevan.index)
    if len(common) == 0:
        raise ValueError('No cell IDs shared between the AnnData object and the '
                         'SCEVAN results.')

    adata.obs['scevan_class'] = pd.NA
    adata.obs.loc[common, 'scevan_class'] = scevan.loc[common, 'class'].values
    print(f'Matched SCEVAN classifications for {len(common):,} of {adata.n_obs:,} cells.')


def classify_clusters(adata):
    """Majority SCEVAN class per cluster, mapped to malignant / non-malignant."""
    mapping = {}
    print('\nCluster-level CNV classification:')

    for cluster, cells in adata.obs.groupby('leiden', observed=True):
        classes = cells['scevan_class'].dropna()
        if classes.empty:
            mapping[str(cluster)] = 'non-malignant'
            print(f'  cluster {cluster:>3}: no SCEVAN predictions -> non-malignant')
            continue

        counts = classes.value_counts()
        scevan_call = counts.idxmax()
        classification = 'malignant' if scevan_call in MALIGNANT_CLASSES else 'non-malignant'
        mapping[str(cluster)] = classification

        print(f'  cluster {cluster:>3}: n={counts.sum():,}; '
              f'SCEVAN={scevan_call} -> {classification}')

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description='CNV-based malignant call from SCEVAN results.')
    parser.add_argument('--input-h5ad', type=Path, required=True)
    parser.add_argument('--scevan-results', type=Path, required=True)
    parser.add_argument('--output-h5ad', type=Path, required=True)
    parser.add_argument('--output-annotations', type=Path, required=True)

    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_h5ad)
    print(f'{adata.n_obs:,} cells x {adata.n_vars:,} genes')
    if 'leiden' not in adata.obs.columns:
        raise ValueError("adata.obs['leiden'] is required.")

    scevan = pd.read_csv(args.scevan_results, index_col=0)
    if 'class' not in scevan.columns:
        raise ValueError("The SCEVAN results file must contain a 'class' column.")
    scevan['class'] = scevan['class'].astype(str)

    add_scevan_predictions(adata, scevan)

    mapping = classify_clusters(adata)
    adata.obs['cnv_annotation'] = pd.Categorical(
        adata.obs['leiden'].astype(str).map(mapping).fillna('non-malignant'),
        categories=['non-malignant', 'malignant'])

    print('\nFinal classification:')
    print(adata.obs['cnv_annotation'].value_counts().to_string())

    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    args.output_annotations.parent.mkdir(parents=True, exist_ok=True)
    adata.write(args.output_h5ad)
    adata.obs[['cnv_annotation']].to_csv(args.output_annotations, sep='\t')

    print(f'\nAnnData:     {args.output_h5ad}')
    print(f'Annotations: {args.output_annotations}')


if __name__ == '__main__':
    main()
