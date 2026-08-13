"""
Per-sample cell type annotation.

Combines the automated consensus annotation from the preprocessing pipeline with
SCEVAN copy-number calls and manual, marker-guided curation of each Leiden
cluster. Non-malignant cells are then reclustered and annotated at finer
resolution.

Run in two stages:

  --stage explore   scores the marker panels, writes violin plots per panel and
                    prints the automatic cluster assignments, so that clusters
                    can be inspected and assigned by eye.
  --stage annotate  applies the manual assignments from the sample config and
                    writes the annotated object.

The manual decisions for each sample live in a small JSON config, so that the
same script produces every sample. Example:

  {
    "entity": "ewing",
    "malignant_clusters": ["8", "10"],
    "cluster_labels": {"10": "Endothelial cells", "12": "Myeloid cells"},
    "non_malignant": {
      "n_pcs": 6,
      "resolution": 1.0,
      "cluster_labels": {"1": "Endothelial cells", "2": "Macrophages"}
    }
  }

The annotation produced here is a starting point. Final cell type and cell state
labels were assigned after cohort-wide integration.

Usage:
  python annotate_sample.py --adata sample.h5ad --config sample.json \
      --scevan scevan_results.csv --output-dir results/ --stage explore
  python annotate_sample.py --adata sample.h5ad --config sample.json \
      --scevan scevan_results.csv --output-dir results/ --stage annotate \
      --output-h5ad sample_annotated.h5ad
"""

import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# Entity-defining markers used to identify malignant cells
MALIGNANT_MARKERS = {
    'ewing': ['RBM11', 'LOXHD1', 'MAPT'],
    'alveolar_rhabdomyosarcoma': ['MYO18B', 'TTN', 'SPATS2L'],
    'embryonal_rhabdomyosarcoma': ['IGF2', 'ZFHX4', 'CNN3'],
    'osteosarcoma': ['RUNX2', 'COL11A1'],
    'clear_cell_sarcoma': ['COL4A2', 'PCP4', 'MDGA2'],
    'desmoplastic_small_round_cell_tumor': ['CACNA2D2'],
    'alveolar_soft_part_sarcoma': ['GPNMB', 'IGFN1', 'PSAP', 'UPP1'],
    'epithelioid_sarcoma': ['STEAP1B', 'MME', 'SORBS2'],
    'synovial_sarcoma': ['TLE1', 'CACNA1C', 'SLIT3'],
    'myoepithelial_carcinoma': ['MUCL1', 'ERBB4'],
    'undifferentiated_sarcoma': ['LAMC1', 'SFRP1'],
}

# Canonical markers used to annotate the non-malignant compartment
NON_MALIGNANT_MARKERS = {
    'Endothelial cells': ['PECAM1', 'CD34', 'VWF', 'KDR'],
    'Macrophages': ['CD68', 'CD163', 'CD14', 'ITGAM'],
    'DC': ['ITGAX', 'CD80', 'CD86'],
    'T cells': ['CD3D', 'CD3G', 'CD8A', 'CD28', 'IL7R'],
    'B cells': ['CD79A', 'RALGPS2', 'MS4A1', 'BANK1', 'IGHM'],
    'NK cells': ['NCAM1', 'FCGR3A', 'KLRK1'],
    'Fibroblasts': ['COL1A1', 'COL3A1', 'DCN', 'FBN1', 'FBLN1'],
    'Pericytes': ['PDGFRB', 'RGS5', 'CD248', 'MCAM'],
    'Immune': ['PTPRC', 'CD4', 'IL7R'],
}


def majority_label_per_cluster(obs, label_col, cluster_col='leiden'):
    """Most frequent value of label_col within each cluster."""
    mapping = {}
    for name, group in obs.groupby(cluster_col, observed=True):
        unique, counts = np.unique(group[label_col].astype(str), return_counts=True)
        mapping[str(name)] = unique[np.argmax(counts)]
    return mapping


def malignant_from_cnv(cluster_to_scevan, malignant_values=('tumor', 'filtered')):
    """
    Clusters called malignant by SCEVAN. Cells SCEVAN could not classify
    ('filtered') sit within tumor clusters and are treated as malignant.
    """
    return [c for c, v in cluster_to_scevan.items() if v in malignant_values]


def present_genes(adata, genes):
    """Marker genes present in the dataset; missing ones are reported."""
    available = set(adata.var_names)
    missing = [g for g in genes if g not in available]
    if missing:
        print(f'   not in dataset: {", ".join(missing)}')
    return [g for g in genes if g in available]


def plot_marker_panels(adata, marker_sets, output_dir, cluster_col='leiden', prefix=''):
    """Violin plot of each marker panel across clusters."""
    os.makedirs(output_dir, exist_ok=True)
    for name, genes in marker_sets.items():
        genes = present_genes(adata, genes)
        if not genes:
            continue
        print(f'-> {name}: {", ".join(genes)}')
        sc.pl.violin(adata, genes, groupby=cluster_col, rotation=70, show=False)
        safe = name.replace(' ', '_').replace('/', '_')
        plt.savefig(os.path.join(output_dir, f'{prefix}{safe}_violin.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()


def recluster(adata, n_pcs, resolution, n_neighbors=15):
    """Recompute PCA, neighborhood graph, UMAP and Leiden on a subset."""
    adata.obs['old_leiden'] = adata.obs['leiden']
    sc.tl.pca(adata, svd_solver='arpack', random_state=42)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                    use_rep='X_pca', random_state=42)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, random_state=42, resolution=resolution)
    print(f'   {adata.obs["leiden"].nunique()} clusters at resolution {resolution}')
    return adata


def main():
    parser = argparse.ArgumentParser(description='Per-sample cell type annotation.')
    parser.add_argument('--adata', required=True,
                        help='Preprocessed .h5ad from the pipeline, with leiden clusters')
    parser.add_argument('--config', required=True,
                        help='JSON with the manual decisions for this sample')
    parser.add_argument('--scevan', default=None,
                        help='SCEVAN results CSV, indexed by cell barcode')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--output-h5ad', default=None,
                        help='Path for the annotated object; required for --stage annotate')
    parser.add_argument('--stage', choices=['explore', 'annotate'], default='explore')
    parser.add_argument('--consensus-col', default='cluster_to_consensus_all',
                        help='Automated consensus annotation from the pipeline')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.config) as f:
        config = json.load(f)

    entity = config['entity']
    if entity not in MALIGNANT_MARKERS:
        raise ValueError(f'Unknown entity: {entity}. '
                         f'Known: {", ".join(sorted(MALIGNANT_MARKERS))}')

    adata = sc.read_h5ad(args.adata)
    print(f'-> {adata.shape[0]} cells x {adata.shape[1]} genes, '
          f'{adata.obs["leiden"].nunique()} clusters')

    # copy-number calls per cluster
    cluster_to_scevan = {}
    if args.scevan is not None:
        scevan = pd.read_csv(args.scevan, index_col=0)
        adata.obs['scevan_res'] = scevan.reindex(adata.obs_names).iloc[:, 0].astype(str)
        cluster_to_scevan = majority_label_per_cluster(adata.obs, 'scevan_res')
        print(f'-> SCEVAN malignant clusters: '
              f'{", ".join(sorted(malignant_from_cnv(cluster_to_scevan)))}')

    if args.stage == 'explore':
        print(f'\n-> Entity markers ({entity})')
        plot_marker_panels(adata, {entity: MALIGNANT_MARKERS[entity]},
                           args.output_dir, prefix='malignant_')
        print('\n-> Non-malignant markers')
        plot_marker_panels(adata, NON_MALIGNANT_MARKERS, args.output_dir)

        print('\n-> Automatic consensus annotation per cluster')
        consensus = majority_label_per_cluster(adata.obs, args.consensus_col)
        for cluster in sorted(consensus, key=lambda x: int(x)):
            cnv = cluster_to_scevan.get(cluster, '-')
            print(f'   cluster {cluster:>3}  consensus: {consensus[cluster]:<20} CNV: {cnv}')
        print(f'\nInspect the plots in {args.output_dir}, then record the cluster '
              f'assignments in {args.config} and rerun with --stage annotate.')
        return

    # --- annotate ---
    if args.output_h5ad is None:
        parser.error('--output-h5ad is required for --stage annotate')

    # start from the pipeline consensus, then apply the manual decisions
    cluster_to_annotation = majority_label_per_cluster(adata.obs, args.consensus_col)

    malignant = set(malignant_from_cnv(cluster_to_scevan))
    malignant |= set(str(c) for c in config.get('malignant_clusters', []))
    for cluster in malignant:
        cluster_to_annotation[cluster] = 'malignant'

    for cluster, label in config.get('cluster_labels', {}).items():
        cluster_to_annotation[str(cluster)] = label

    adata.obs['final_annotation'] = adata.obs['leiden'].astype(str).map(cluster_to_annotation)
    print('-> First-pass annotation')
    print(adata.obs['final_annotation'].value_counts().to_string())

    # recluster and annotate the non-malignant compartment
    nm_config = config.get('non_malignant')
    if nm_config:
        print('\n-> Reclustering non-malignant cells')
        adata_nm = adata[adata.obs['final_annotation'] != 'malignant'].copy()
        adata_nm = recluster(adata_nm, nm_config['n_pcs'], nm_config['resolution'])

        if args.stage == 'annotate' and not nm_config.get('cluster_labels'):
            plot_marker_panels(adata_nm, NON_MALIGNANT_MARKERS,
                               args.output_dir, prefix='non_malignant_')
            print(f'   no non-malignant cluster labels in the config yet; '
                  f'violin plots written to {args.output_dir}')
        else:
            labels = {str(k): v for k, v in nm_config['cluster_labels'].items()}
            adata_nm.obs['final_annotation'] = (
                adata_nm.obs['leiden'].astype(str).map(labels).fillna('unclear'))
            adata.obs['final_annotation'] = adata.obs['final_annotation'].astype(str)
            adata.obs.loc[adata_nm.obs_names, 'final_annotation'] = \
                adata_nm.obs['final_annotation']
            print('\n-> Final annotation')
            print(adata.obs['final_annotation'].value_counts().to_string())

    adata.write(args.output_h5ad)
    adata.obs['final_annotation'].to_csv(
        os.path.join(args.output_dir, 'final_annotation.csv'), sep='\t')
    print(f'\nDone. Annotated object written to {args.output_h5ad}')


if __name__ == '__main__':
    main()
