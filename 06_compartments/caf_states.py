"""
CAF state characterization.

Generates the evidence used to annotate cancer-associated fibroblast states:
subclustering of the stromal compartment, marker panel scores per cluster,
differential expression, and marker expression plots. The cell states themselves
were assigned by manual review of these outputs.

Marker panels are those of Cords et al. (Nat Commun 2023, 14:4294), in a short
form of one to three defining genes per state and a long form used for the
expression matrixplot.

Input: .h5ad of the stromal compartment, log-normalized with the full gene set
in .raw and an existing UMAP embedding.

Output: figures and the per-cluster score table in --output-dir, and the
reclustered object at --output-h5ad.

Usage:
  python caf_states.py --adata stromal.h5ad --output-dir figures/ \
      --output-h5ad stromal_reclustered.h5ad --resolution 1.0
"""

import argparse
import os
import harmonypy as hm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

# Defining markers per CAF state (Cords et al. 2023)
MARKERS_SHORT = {
    'Pericyte': ['RGS5'],
    'vCAFs': ['MCAM', 'MYH11'],
    'tCAFs': ['MME', 'GAPDH'],
    'hsp_tCAFs': ['HSPH1'],
    'ifnCAFs': ['IDO1'],
    'mCAFs': ['MMP11', 'POSTN', 'COL1A1'],
    'iCAFs': ['PLA2G2A', 'CFD', 'CD34'],
    'dCAFs': ['MKI67', 'TUBA1B'],
    'rCAFs': ['CCL21', 'CCL19'],
    'apCAFs': ['CD74', 'HLA-DRA'],
}

# Extended panels used for the expression matrixplot
MARKERS_LONG = {
    'Pericyte': ['RGS5'],
    'vCAFs': ['MCAM', 'NOTCH3', 'COL18A1', 'MYH11'],
    'tCAFs': ['MME', 'NDRG1', 'ENO1', 'GAPDH', 'VEGFA'],
    'ifnCAFs': ['CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'IL32'],
    'hsp_tCAFs': ['HSPH1', 'HSP90AA1'],
    'mCAFs': ['MMP11', 'COL1A1', 'POSTN', 'COL10A1', 'COL11A1', 'COL8A1',
              'COL1A2', 'COL12A1', 'COL3A1', 'COL5A2', 'COMP', 'LRRC15',
              'LRRC17', 'ASPN', 'SULF1', 'VCAN', 'INHBA', 'MGP', 'BGN'],
    'iCAFs': ['PLA2G2A', 'CFD', 'CD34', 'C3', 'CXCL12', 'CXCL14', 'IL6'],
    'dCAFs': ['MKI67', 'TUBA1B'],
    'rCAFs': ['CCL21', 'CCL19'],
    'apCAFs': ['HLA-DRA', 'HLA-DRB1', 'CD74'],
}

def integrate(adata, n_pcs=50, resolution=0.5, n_neighbors=15, cluster_key="leiden"):
    """PCA, Harmony integration, neighbourhood graph, UMAP and Leiden."""
    sc.tl.pca(adata, svd_solver='arpack', random_state=42)
    adata.obs['10x_chip'] = adata.obs['10x_chip'].astype(str)
    harmony_out = hm.run_harmony(adata.obsm["X_pca"].copy(), adata.obs, ['sample', '10x_chip'])
    adata.obsm["X_pca_harmony"] = harmony_out.Z_corr
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                    use_rep='X_pca_harmony', random_state=42)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, random_state=42, resolution=resolution, key_added=cluster_key)
    print(f'-> {adata.obs[cluster_key].nunique()} clusters at resolution {resolution}')
    return adata


def score_marker_panels(adata, marker_sets, cluster_key, output_dir):
    """
    Score each marker panel per cell, then take the mean per cluster. Each
    cluster is labelled with its highest-scoring panel as a starting point for
    manual review.
    """
    available = adata.var_names
    scores = []

    for state, genes in marker_sets.items():
        present = [g for g in genes if g in available]
        if not present:
            print(f'   no genes present for {state}, skipping')
            continue
        sc.tl.score_genes(adata, present, score_name=state)
        scores.append(state)

    cluster_means = adata.obs.groupby(cluster_key)[scores].mean()
    adata.obs['top_scoring_state'] = adata.obs[cluster_key].map(
        cluster_means.idxmax(axis=1))

    cluster_means.to_csv(os.path.join(output_dir, 'cluster_marker_scores.csv'))

    sc.pl.umap(adata, color=scores, legend_fontsize=10, legend_fontoutline=3,
               wspace=0.2, cmap='bwr', vmin=-4, show=False)
    plt.savefig(os.path.join(output_dir, 'marker_scores_umap.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10, 8))
    sns.heatmap(cluster_means.T, cmap='bwr', vmin=-1.5, vmax=1.5,
                linewidths=.5, xticklabels=True)
    plt.title('Mean marker panel score per cluster')
    plt.savefig(os.path.join(output_dir, 'marker_scores_heatmap.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    return cluster_means


def marker_expression(adata, marker_sets, cluster_key, output_dir):
    """Matrixplot of scaled marker expression per cluster."""
    #adata_raw = adata.raw.to_adata()
    adata_raw = adata
    adata_raw.obs[cluster_key] = adata.obs[cluster_key]
    adata_raw.obsm['X_umap'] = adata.obsm['X_umap']
    adata_raw.layers['scaled'] = sc.pp.scale(adata_raw, copy=True).X
    sc.tl.dendrogram(adata_raw, groupby=cluster_key)

    present = {state: [g for g in genes if g in adata_raw.var_names]
               for state, genes in marker_sets.items()}
    present = {k: v for k, v in present.items() if v}

    sc.pl.matrixplot(adata_raw, present, cluster_key, dendrogram=True,
                     colorbar_title='mean z-score', layer='scaled',
                     vmin=-2, vmax=2, cmap='RdBu_r', swap_axes=True, show=False)
    plt.savefig(os.path.join(output_dir, 'marker_expression_matrixplot.png'),
                dpi=300, bbox_inches='tight')
    plt.close()


def differential_expression(adata, cluster_key, output_dir):
    """Rank marker genes per cluster and plot the top genes."""
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, swap_axes=False, show=False)
    plt.savefig(os.path.join(output_dir, 'DEG_dotplot.png'), dpi=300, bbox_inches='tight')
    plt.close()

    sc.get.rank_genes_groups_df(adata, group=None).to_csv(
        os.path.join(output_dir, 'DEG_per_cluster.csv'), index=False)


def main():
    parser = argparse.ArgumentParser(description='CAF state characterization.')
    parser.add_argument('--adata', required=True,
                        help='.h5ad of the stromal compartment')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--output-h5ad', required=True)
    parser.add_argument('--resolution', type=float, default=1.0)
    parser.add_argument('--cluster-key', default='leiden')
    parser.add_argument('--n-pcs', type=int, default=50)

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir
    sc.settings.set_figure_params(dpi=80, dpi_save=300, facecolor='white')

    adata = sc.read_h5ad(args.adata)
    print(f'{adata.shape[0]} cells x {adata.shape[1]} genes')

    print('-> Reintegrating and clustering')
    adata = integrate(adata, n_pcs=args.n_pcs, resolution=args.resolution, cluster_key=args.cluster_key)
    cluster_key = args.cluster_key

    sc.pl.umap(adata, color=[cluster_key], legend_loc='on data',
               legend_fontsize=12, legend_fontoutline=4, show=False)
    plt.savefig(os.path.join(args.output_dir, 'umap_clusters.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    print('-> Differential expression per cluster')
    differential_expression(adata, cluster_key, args.output_dir)

    print('-> Scoring CAF marker panels')
    score_marker_panels(adata, MARKERS_SHORT, cluster_key, args.output_dir)

    print('-> Marker expression per cluster')
    marker_expression(adata, MARKERS_LONG, cluster_key, args.output_dir)

    adata.write(args.output_h5ad)
    print(f'\nDone. Figures in {args.output_dir}, object at {args.output_h5ad}')


if __name__ == '__main__':
    main()
