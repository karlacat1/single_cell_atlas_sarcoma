"""
T cell state characterization.

Generates the evidence used to annotate CD4+ and CD8+ T cell states: Harmony
reintegration and clustering, functional signature scores, scores for the
reference cell state panels, differential expression and marker expression
plots. The cell states themselves were assigned by manual review of these
outputs.

Reference cell state panels are the top marker genes of the pan-cancer T cell
atlas of Chu et al. (Nat Med 2023, 29:1550); all panels are scored so that
states absent from this cohort are visible as such.

Input: .h5ad of the CD4+ or CD8+ compartment, log-normalized, with obs columns
sample and 10x_chip.

Output: figures and the per-cluster score table in --output-dir, and the
reintegrated object at --output-h5ad.

Usage:
  python tcell_states.py --adata cd8.h5ad --compartment CD8 \
      --output-dir figures/cd8/ --output-h5ad cd8_reclustered.h5ad
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
import scanpy.external as sce
import seaborn as sns

# Functional signatures, scored per cell and averaged per cluster
SIGNATURES = {
    'CD4': {
        'Naive': ['IL7R', 'CCR7', 'SELL', 'FOXP1', 'KLF2', 'KLF3', 'LEF1', 'TCF7',
                  'ACTN1', 'BTG1', 'BTG2', 'TOB1'],
        'Activation/Effector': ['FAS', 'CD44', 'CD69', 'CD38', 'NKG7', 'KLRB1', 'KLRD1',
                                'KLRG1', 'CX3CR1', 'CD300A', 'FGFBP2', 'ID2', 'ID3',
                                'PRDM1', 'RUNX3', 'TBX21', 'ZEB2', 'BATF', 'NR4A1',
                                'NR4A2', 'HOPX', 'FOS', 'FOSB', 'FOSL2', 'JUN', 'JUNB',
                                'JUND', 'STAT1', 'STAT3', 'EOMES', 'AHR'],
        'Exhaustion': ['PDCD1', 'LAYN', 'HAVCR2', 'LAG3', 'CTLA4', 'TIGIT', 'TOX',
                       'VSIR', 'BTLA', 'ENTPD1'],
    },
    'CD8': {
        'Naive': ['IL7R', 'CCR7', 'SELL', 'FOXO1', 'KLF2', 'KLF3', 'LEF1', 'TCF7',
                  'ACTN1', 'FOXP1'],
        'Activation/Effector': ['FAS', 'FASLG', 'CD44', 'CD69', 'CD38', 'NKG7', 'KLRB1',
                                'KLRD1', 'KLRF1', 'KLRG1', 'KLRK1', 'FCGR3A', 'CX3CR1',
                                'CD300A', 'FGFBP2', 'ID2', 'ID3', 'PRDM1', 'RUNX3',
                                'TBX21', 'ZEB2', 'BATF', 'IRF4', 'NR4A1', 'NR4A2',
                                'NR4A3', 'PBX3', 'ZNF683', 'HOPX', 'FOS', 'FOSB', 'JUN',
                                'JUNB', 'JUND', 'STAT1', 'STAT2', 'STAT5A', 'STAT6',
                                'STAT4', 'EOMES'],
        'Cytotoxicity': ['GZMA', 'GZMB', 'GZMH', 'GZMK', 'GNLY', 'PRF1', 'IFNG', 'TNF',
                         'SERPINB1', 'SERPINB6', 'SERPINB9', 'CTSA', 'CTSB', 'CTSC',
                         'CTSD', 'CTSW', 'CST3', 'CST7', 'CSTB', 'LAMP1', 'LAMP3',
                         'CAPN2'],
        'Exhaustion': ['PDCD1', 'LAYN', 'HAVCR2', 'LAG3', 'CD244', 'CTLA4', 'LILRB1',
                       'TIGIT', 'TOX', 'VSIR', 'BTLA', 'ENTPD1', 'CD160', 'LAIR1'],
    },
}

# Reference cell state panels from the pan-cancer T cell atlas (Chu et al. 2023)
REFERENCE_STATES = {
    'CD4': {
        'CD4_c0_Tcm': ['IL7R', 'GPR183', 'CD69'],
        'CD4_c1_Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF4'],
        'CD4_c2_Tn': ['RPL31', 'RPL21'],
        'CD4_c3_Tfh': ['CXCL13', 'PDCD1', 'TOX', 'ICOS', 'BCL6'],
        'CD4_c4_Tstr': ['NR4A1', 'BAG3', 'FOS', 'JUN'],
        'CD4_c5_CTL': ['IFNG', 'GZMA', 'GZMH', 'GZMB', 'GZMK', 'NKG7', 'PRF1'],
        'CD4_c6_Tn_FHIT': ['FHIT', 'CCR7', 'LEF1', 'SELL', 'TCF7'],
        'CD4_c7_Tn_TCEA3': ['CCR7', 'LEF1', 'SELL', 'TCF7', 'TCEA3', 'CDC25B'],
        'CD4_c8_Th17': ['IL17F', 'IL17A', 'RORA', 'KLRB1', 'CCR6'],
        'CD4_c9_Tn_TCF7_SLC40A1': ['SLC40A1'],
        'CD4_c10_Tn_LEF1_ANKRD55': ['ANKRD55'],
        'CD4_c11_Tisg': ['ISG15', 'IFI44L', 'IFIT1'],
    },
    'CD8': {
        'CD8_c0_t-Teff': ['GZMK', 'GZMB', 'PRF1', 'CD44', 'CD69'],
        'CD8_c1_Tex': ['FAS', 'FASLG', 'PDCD1', 'LAG3', 'CTLA4', 'TOX', 'TIGIT'],
        'CD8_c2_Teff': ['FGFBP2', 'GZMH', 'GNLY'],
        'CD8_c3_Tn': ['CCR7', 'SELL'],
        'CD8_c4_Tstr': ['NR4A1', 'BAG3', 'HSPA1A', 'HSPA1B'],
        'CD8_c5_Tisg': ['IFIT1', 'MX1'],
        'CD8_c6_Tcm': ['DKK3', 'CCR4', 'EOMES'],
        'CD8_c7_p-Tex': ['EOMES', 'CNN2', 'LIMD2', 'CD27'],
        'CD8_c8_Teff_KLRG1': ['LGR6'],
        'CD8_c9_Tsen': ['KLRC4'],
        'CD8_c10_Teff_CD244': ['CD244'],
        'CD8_c11_Teff_SEMA4A': ['SEMA4A'],
        'CD8_c12_Trm': ['ITGA1', 'KLRB1', 'PRDM1'],
        'CD8_c13_Tn_TCF7': ['CCR7', 'SELL', 'TCF7'],
    },
}


def integrate(adata, n_pcs=50, resolution=0.5, n_neighbors=15):
    """PCA, Harmony integration, neighbourhood graph, UMAP and Leiden."""
    sc.tl.pca(adata, svd_solver='arpack', random_state=42)
    adata.obs['10x_chip'] = adata.obs['10x_chip'].astype(str)
    sce.pp.harmony_integrate(adata, ['sample', '10x_chip'])
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                    use_rep='X_pca_harmony', random_state=42)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, random_state=42, resolution=resolution)
    print(f'-> {adata.obs["leiden"].nunique()} clusters at resolution {resolution}')
    return adata


def score_panels(adata, panels, cluster_key, output_dir, prefix, label_clusters=False):
    """
    Score each panel per cell and average per cluster. Panels with no genes
    present in the dataset are reported and skipped.
    """
    scores = []
    for name, genes in panels.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            print(f'   no genes present for {name}, skipping')
            continue
        sc.tl.score_genes(adata, present, score_name=name)
        scores.append(name)

    cluster_means = adata.obs.groupby(cluster_key)[scores].mean()
    cluster_means.to_csv(os.path.join(output_dir, f'{prefix}_cluster_scores.csv'))

    if label_clusters:
        adata.obs[f'{prefix}_top_scoring'] = adata.obs[cluster_key].map(
            cluster_means.idxmax(axis=1))

    sc.pl.umap(adata, color=scores, legend_fontsize=10, legend_fontoutline=3,
               wspace=0.3, cmap='bwr', ncols=4, show=False)
    plt.savefig(os.path.join(output_dir, f'{prefix}_scores_umap.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    sc.pl.violin(adata, keys=scores, groupby=cluster_key, show=False)
    plt.savefig(os.path.join(output_dir, f'{prefix}_scores_violin.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10, 8))
    sns.heatmap(cluster_means.T, cmap='bwr', center=0, linewidths=.5, xticklabels=True)
    plt.title(f'Mean {prefix} score per cluster')
    plt.savefig(os.path.join(output_dir, f'{prefix}_scores_heatmap.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    return cluster_means


def differential_expression(adata, cluster_key, output_dir):
    """Rank marker genes per cluster and plot the top genes."""
    sc.tl.dendrogram(adata, groupby=cluster_key)
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
    plt.savefig(os.path.join(output_dir, 'DEG_dotplot.png'), dpi=300, bbox_inches='tight')
    plt.close()

    sc.get.rank_genes_groups_df(adata, group=None).to_csv(
        os.path.join(output_dir, 'DEG_per_cluster.csv'), index=False)


def marker_expression(adata, panels, cluster_key, output_dir):
    """Matrixplot of reference panel expression per cluster."""
    present = {name: [g for g in genes if g in adata.var_names]
               for name, genes in panels.items()}
    present = {k: v for k, v in present.items() if v}

    sc.pl.matrixplot(adata, present, cluster_key, dendrogram=True, log=True,
                     cmap='RdBu_r', vmin=-1, vmax=1, swap_axes=False, show=False)
    plt.savefig(os.path.join(output_dir, 'reference_markers_matrixplot.png'),
                dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='T cell state characterization.')
    parser.add_argument('--adata', required=True,
                        help='.h5ad of the CD4+ or CD8+ compartment')
    parser.add_argument('--compartment', required=True, choices=['CD4', 'CD8'])
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--output-h5ad', required=True)
    parser.add_argument('--resolution', type=float, default=0.5)
    parser.add_argument('--n-pcs', type=int, default=50)
    parser.add_argument('--cluster-key', default='leiden')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir
    sc.settings.set_figure_params(dpi=80, dpi_save=300, facecolor='white')

    adata = sc.read_h5ad(args.adata)
    print(f'{adata.shape[0]} cells x {adata.shape[1]} genes')

    print('-> Reintegrating and clustering')
    adata = integrate(adata, n_pcs=args.n_pcs, resolution=args.resolution)

    sc.pl.umap(adata, color=[args.cluster_key], legend_loc='on data',
               legend_fontsize=10, legend_fontoutline=4, show=False)
    plt.savefig(os.path.join(args.output_dir, 'umap_clusters.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    print('-> Differential expression per cluster')
    differential_expression(adata, args.cluster_key, args.output_dir)

    print('-> Scoring functional signatures')
    score_panels(adata, SIGNATURES[args.compartment], args.cluster_key,
                 args.output_dir, 'signature')

    print('-> Scoring reference cell state panels')
    score_panels(adata, REFERENCE_STATES[args.compartment], args.cluster_key,
                 args.output_dir, 'reference_state', label_clusters=True)

    print('-> Reference marker expression per cluster')
    marker_expression(adata, REFERENCE_STATES[args.compartment],
                      args.cluster_key, args.output_dir)

    adata.write(args.output_h5ad)
    print(f'\nDone. Figures in {args.output_dir}, object at {args.output_h5ad}')


if __name__ == '__main__':
    main()
