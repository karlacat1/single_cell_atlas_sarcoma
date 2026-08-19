"""
Macrophage state characterization.

Generates the evidence used to annotate macrophage cell states: differential
expression per cluster, reference marker panel scores, GO term enrichment and
Reactome pathway enrichment. The cell states themselves were assigned by manual
review of these outputs.

Ribosomal genes (RPS, RPL) are removed and the clustering recomputed before
scoring, since ribosomal content otherwise dominates the macrophage subclusters.

Input: .h5ad of the macrophage compartment, log-normalized with the full gene
set in .raw and obs columns sample and 10x_chip.

Output: figures and enrichment tables in --output-dir, and the reclustered
object at --output-h5ad.

Usage:
  python macrophage_states.py --adata macrophages.h5ad --output-dir figures/ \
      --output-h5ad macrophages_noribo.h5ad \
      --marker-table markers_guimaraes.xlsx --pathway-table pathways.xlsx
"""

import argparse
import os
import harmonypy as hm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import gseapy

# Macrophage-relevant GO terms used to filter the enrichment results
MACROPHAGE_GO_TERMS = [
    'GO:0042116',  # macrophage activation
    'GO:0043030',  # regulation of macrophage activation
    'GO:0002472',  # macrophage antigen processing and presentation
]

# TAM markers used as an orthogonal check on the suppressive states
TAM_MARKERS = ['CD163', 'MRC1', 'ARG1', 'VEGFA', 'IL10', 'TGFB1', 'MMP9',
               'MMP2', 'CCL2', 'CSF1R', 'CD274', 'SIRPA', 'SPP1', 'MARCO']


def remove_ribosomal_genes(adata):
    """Drop ribosomal genes and report how many were removed."""
    ribo = adata.var_names.str.startswith(('RPS', 'RPL'))
    print(f'-> Removing {ribo.sum()} ribosomal genes of {adata.n_vars}')
    return adata[:, ~ribo].copy()


def recluster(adata, n_pcs=50, resolution=1.0, n_neighbors=15):
    """PCA, Harmony integration, neighbourhood graph, UMAP and Leiden."""
    sc.tl.pca(adata, svd_solver='arpack', random_state=42)
    adata.obs['10x_chip'] = adata.obs['10x_chip'].astype(str)
    harmony_out = hm.run_harmony(adata.obsm["X_pca"].copy(), adata.obs, ['sample', '10x_chip'])
    adata.obsm["X_pca_harmony"] = harmony_out.Z_corr
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                    use_rep='X_pca_harmony', random_state=42)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, resolution=resolution, random_state=42)
    print(f'-> {adata.obs["leiden"].nunique()} clusters')
    return adata


def differential_expression(adata, output_dir):
    """Rank marker genes per cluster and plot the top genes."""
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
    sc.tl.dendrogram(adata, groupby='leiden')

    sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden', n_genes=5, show=False)
    plt.savefig(os.path.join(output_dir, 'DEG_dotplot.png'), dpi=300, bbox_inches='tight')
    plt.close()

    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, figsize=(10, 7), show=False)
    plt.savefig(os.path.join(output_dir, 'DEG_per_cluster.png'), dpi=300, bbox_inches='tight')
    plt.close()


def score_marker_panels(adata, marker_sets, output_dir, prefix):
    """Score each marker panel and plot the scores on the UMAP."""
    for name, genes in marker_sets.items():
        genes = [g for g in genes if g in adata.var_names]
        if not genes:
            print(f'   no genes present for {name}, skipping')
            continue
        sc.tl.score_genes(adata, gene_list=genes, score_name=f'{name}_score')

    scores = [f'{n}_score' for n in marker_sets if f'{n}_score' in adata.obs.columns]
    if scores:
        sc.pl.umap(adata, color=['leiden'] + scores, ncols=5, wspace=0.2, show=False)
        plt.savefig(os.path.join(output_dir, f'{prefix}_scores.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()


def enrichment_per_cluster(adata, gene_sets, output_dir, prefix, top_n=200):
    """Enrichment of each cluster's upregulated genes against a gene set library."""
    results = []

    for cluster in sorted(set(adata.obs['leiden']), key=int):
        degs = sc.get.rank_genes_groups_df(adata, group=cluster)
        degs = degs[(degs['pvals_adj'] < 0.05) & (degs['logfoldchanges'] > 0)]
        genes = degs.sort_values('logfoldchanges', ascending=False)['names'].head(top_n).tolist()
        if not genes:
            continue

        enr = gseapy.enrichr(gene_list=genes, gene_sets=gene_sets, outdir=None)
        res = enr.results
        res['cluster'] = cluster
        results.append(res)

    if not results:
        return pd.DataFrame()

    all_res = pd.concat(results, ignore_index=True)
    all_res.to_csv(os.path.join(output_dir, f'{prefix}_enrichment.tsv'),
                   sep='\t', index=False)
    return all_res


def enrichment_heatmap(enrichment, terms, output_path, title):
    """Heatmap of adjusted p-values for selected terms across clusters."""
    subset = enrichment[enrichment['Term'].isin(terms)]
    if subset.empty:
        print(f'   no matching terms for {title}, skipping heatmap')
        return

    heatmap_data = subset.pivot_table(index='Term', columns='cluster',
                                      values='Adjusted P-value', aggfunc='first')

    g = sns.clustermap(-np.log10(heatmap_data.fillna(1)), cmap='coolwarm',
                       linewidths=.5, figsize=(10, max(4, len(heatmap_data) * 0.3)),
                       cbar_kws={'label': '-log10(adjusted P)'})
    g.figure.suptitle(title)
    g.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Macrophage state characterization.')
    parser.add_argument('--adata', required=True,
                        help='.h5ad of the macrophage compartment')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--output-h5ad', required=True)
    parser.add_argument('--marker-table', default=None,
                        help='Excel table of reference marker panels, with columns '
                             'cluster and gene')
    parser.add_argument('--pathway-table', default=None,
                        help='Excel table of reference pathways, with a Description column')
    parser.add_argument('--go-obo', default=None,
                        help='GO basic OBO file, used to resolve GO term names')
    parser.add_argument('--resolution', type=float, default=1.0)
    parser.add_argument('--n-pcs', type=int, default=50)

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    sc.settings.figdir = args.output_dir
    sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300,
                                  figsize=(5, 4), format='png', frameon=False)

    adata = sc.read_h5ad(args.adata)
    print(f'{adata.shape[0]} cells x {adata.shape[1]} genes')

    print('-> Removing ribosomal genes and reclustering')
    adata = remove_ribosomal_genes(adata)
    adata = recluster(adata, n_pcs=args.n_pcs, resolution=args.resolution)

    sc.pl.umap(adata, color='leiden', show=False)
    plt.savefig(os.path.join(args.output_dir, 'umap_leiden.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    print('-> Differential expression per cluster')
    differential_expression(adata, args.output_dir)

    # reference marker panels
    if args.marker_table:
        print('-> Scoring reference marker panels')
        markers = pd.read_excel(args.marker_table)
        marker_sets = markers.groupby('cluster')['gene'].apply(list).to_dict()
        score_marker_panels(adata, marker_sets, args.output_dir, 'reference_markers')

    print('-> Scoring TAM markers')
    score_marker_panels(adata, {'TAM': TAM_MARKERS}, args.output_dir, 'TAM')

    # GO term enrichment
    print('-> GO term enrichment')
    go_enrichment = enrichment_per_cluster(
        adata, ['GO_Biological_Process_2023'], args.output_dir, 'GO')
    if not go_enrichment.empty:
        macrophage_terms = go_enrichment[
            go_enrichment['Term'].str.contains('macrophage', case=False)]['Term'].unique()
        enrichment_heatmap(go_enrichment, macrophage_terms,
                           os.path.join(args.output_dir, 'GO_macrophage_heatmap.png'),
                           'Macrophage-related GO terms')

    # Reactome pathway enrichment
    print('-> Reactome pathway enrichment')
    pathway_enrichment = enrichment_per_cluster(
        adata, ['Reactome_Pathways_2024'], args.output_dir, 'pathways')
    if not pathway_enrichment.empty and args.pathway_table:
        reference_terms = pd.read_excel(args.pathway_table)['Description'].tolist()
        enrichment_heatmap(pathway_enrichment, reference_terms,
                           os.path.join(args.output_dir, 'pathway_heatmap.png'),
                           'Reference pathways')

    adata.write(args.output_h5ad)
    print(f'\nDone. Figures in {args.output_dir}, object at {args.output_h5ad}')


if __name__ == '__main__':
    main()
