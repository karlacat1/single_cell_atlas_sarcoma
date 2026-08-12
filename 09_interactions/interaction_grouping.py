"""
Group samples by their cell-cell communication landscape.

Each sample is represented by a vector over (sender cell type, receiver cell type,
ligand complex, receptor complex) features, scored by how strongly and how
specifically each interaction is expressed. Samples are then embedded, clustered,
and the interactions driving each cluster are identified.

Steps:
  1. Load per-sample LIANA results and keep interactions passing the magnitude
     and specificity rank thresholds.
  2. Score each retained interaction as (1 - magnitude_rank) * (1 - specificity_rank),
     so that interactions which are both strongly and specifically expressed
     receive higher values.
  3. Assemble a sample-by-interaction matrix, with absent interactions set to zero
     and interactions with no variance across samples removed.
  4. Scale without mean-centering (the matrix is sparse) and reduce to
     n_components by truncated SVD, retaining the first n_pcs.
  5. Build a nearest-neighbour graph on those components and cluster the samples
     with Leiden.
  6. Test each interaction for differential activity across clusters
     (Kruskal-Wallis), then one-versus-rest per cluster (Mann-Whitney U), with
     Benjamini-Hochberg correction.

Input:
  An .h5ad file whose .obs contains a sample column, plus per-sample LIANA result
  tables. LIANA results are read either from adata.uns (keys of the form
  '<uns_prefix><sample>') or from a directory of pickled DataFrames
  (see --liana-dir).

  Each LIANA table must contain the columns: source, target, ligand_complex,
  receptor_complex, magnitude_rank, specificity_rank.

Output (written to --output-dir):
  interaction_matrix.csv        sample-by-interaction score matrix
  sample_embedding.csv          SVD and UMAP coordinates with cluster assignment
  kruskal_results.csv           across-cluster test per interaction
  one_vs_rest_results.csv       per-cluster one-versus-rest test per interaction
  top_interactions_heatmap.csv  mean score per cluster for the top interactions
  top_interactions_heatmap.png  z-scored heatmap of those interactions

Usage:
  python interaction_grouping.py --adata atlas.h5ad --output-dir results/
  python interaction_grouping.py --adata atlas.h5ad --liana-dir ccc/per_sample/ \
      --group cell_states_detailed --output-dir results/
"""

import argparse
import os
import pickle
import re

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import kruskal, mannwhitneyu
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

REQUIRED_COLUMNS = ['source', 'target', 'ligand_complex', 'receptor_complex',
                    'magnitude_rank', 'specificity_rank']


def load_liana_results(adata, sample_col, uns_prefix, liana_dir, group):
    """
    Collect the per-sample LIANA result tables.

    Results are taken from adata.uns when a key '<uns_prefix><sample>' exists,
    and otherwise from <liana_dir>/<sample>/uns_liana_<sample>_<group>.pkl.

    :return: dict mapping sample name to LIANA DataFrame
    """
    samples = adata.obs[sample_col].astype(str).unique()
    results = {}

    for sample in samples:
        uns_key = f'{uns_prefix}{sample}'

        if uns_key in adata.uns:
            results[sample] = pd.DataFrame(adata.uns[uns_key]).copy()
            continue

        if liana_dir is not None:
            path = os.path.join(liana_dir, sample, f'uns_liana_{sample}_{group}.pkl')
            if os.path.exists(path):
                with open(path, 'rb') as f:
                    results[sample] = pickle.load(f)
                results[sample] = pd.DataFrame(results[sample]).copy()
                continue

        print(f'No cell-cell communication results for {sample}')

    if not results:
        raise ValueError('No LIANA results found. Check --uns-prefix and --liana-dir.')

    for sample, df in results.items():
        missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(f'LIANA results for {sample} are missing columns: {missing}')

    print(f'Loaded LIANA results for {len(results)} of {len(samples)} samples')
    return results


def build_interaction_matrix(liana_results, threshold):
    """
    Build the sample-by-interaction score matrix.

    Interactions are retained when both the magnitude and specificity ranks are
    at or below `threshold`, and scored as the product of (1 - magnitude_rank)
    and (1 - specificity_rank).

    :return: DataFrame, rows = samples, columns = interaction features
    """
    all_rows = []

    for sample, liana_res in liana_results.items():
        liana_res = liana_res.copy()

        # feature identity: sender__receiver|ligand__receptor
        liana_res['feature'] = (liana_res['source'] + '__' + liana_res['target'] + '|' +
                                liana_res['ligand_complex'] + '__' + liana_res['receptor_complex'])

        # keep only interactions passing both rank thresholds
        keep = ((liana_res['magnitude_rank'] <= threshold) &
                (liana_res['specificity_rank'] <= threshold))
        liana_res = liana_res[keep]

        # combined score: high when the interaction is both strong and specific
        liana_res['score'] = ((1 - liana_res['magnitude_rank']) *
                              (1 - liana_res['specificity_rank']))
        liana_res['sample'] = sample

        all_rows.append(liana_res[['sample', 'feature', 'score']])

    long_df = pd.concat(all_rows, ignore_index=True)

    interaction_matrix = long_df.pivot_table(
        index='sample',
        columns='feature',
        values='score',
        aggfunc='max',      # in case of duplicate features within a sample
        fill_value=0
    )

    # drop interactions with no variance across samples
    nonzero_var = interaction_matrix.columns[interaction_matrix.var(axis=0) > 0]
    interaction_matrix = interaction_matrix[nonzero_var]

    print(f'Interaction matrix: {interaction_matrix.shape[0]} samples x '
          f'{interaction_matrix.shape[1]} interactions')
    return interaction_matrix


def reduce_dimensions(interaction_matrix, n_components, random_state=42):
    """
    Scale without mean-centering and reduce with truncated SVD.

    :return: (DataFrame of components, fitted TruncatedSVD)
    """
    # with_mean=False because the matrix is sparse (mostly zeros)
    scaler = StandardScaler(with_mean=False)
    x_scaled = scaler.fit_transform(interaction_matrix.values)

    n_components = min(n_components, min(interaction_matrix.shape) - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    coords = svd.fit_transform(x_scaled)

    components = pd.DataFrame(
        coords,
        index=interaction_matrix.index,
        columns=[f'PC{i + 1}' for i in range(n_components)]
    )
    print(f'Truncated SVD: {n_components} components, '
          f'{svd.explained_variance_ratio_.sum():.1%} of variance explained')
    return components, svd


def cluster_samples(components, n_pcs, n_neighbors, resolution, random_state=42):
    """
    Build a nearest-neighbour graph on the first n_pcs components and cluster
    the samples with Leiden.

    :return: Series of cluster labels indexed by sample
    """
    n_pcs = min(n_pcs, components.shape[1])

    # wrap the coordinates in a temporary AnnData to use scanpy's graph and Leiden
    tmp = ad.AnnData(X=components.iloc[:, :n_pcs].values)
    tmp.obs.index = components.index

    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, use_rep='X', random_state=random_state)
    sc.tl.leiden(tmp, resolution=resolution, random_state=random_state)

    clusters = tmp.obs['leiden'].rename('cluster')
    clusters.index = components.index

    print(f'Leiden clustering: {clusters.nunique()} interaction groups')
    print(clusters.value_counts().sort_index().to_string())
    return clusters


def compute_umap(components, n_pcs, n_neighbors, min_dist=0.1, random_state=42):
    """
    UMAP embedding of the samples in interaction space.

    :return: DataFrame with UMAP1 and UMAP2 indexed by sample
    """
    import umap

    n_pcs = min(n_pcs, components.shape[1])
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist,
                        n_components=2, random_state=random_state)
    coords = reducer.fit_transform(components.iloc[:, :n_pcs].values)

    return pd.DataFrame(coords, index=components.index, columns=['UMAP1', 'UMAP2'])


def test_across_clusters(interaction_matrix, clusters):
    """
    Kruskal-Wallis test per interaction across all clusters, with
    Benjamini-Hochberg correction.

    Interactions are skipped when any cluster has fewer than two samples, when
    all values are identical, or when any cluster has zero variance.
    """
    unique_clusters = sorted(clusters.unique())
    results = []

    for feature in interaction_matrix.columns:
        groups = [interaction_matrix.loc[clusters == c, feature].values
                  for c in unique_clusters]

        if any(len(g) < 2 for g in groups):
            continue
        if len(np.unique(np.concatenate(groups))) < 2:
            continue
        if any(len(np.unique(g)) < 2 for g in groups):
            continue

        try:
            h, p = kruskal(*groups)
        except Exception:
            continue
        if np.isnan(h) or np.isnan(p):
            continue
        results.append({'feature': feature, 'H': h, 'pval': p})

    results_df = pd.DataFrame(results).dropna(subset=['pval'])
    if len(results_df):
        results_df['pval_adj'] = multipletests(results_df['pval'], method='fdr_bh')[1]
        results_df = results_df.sort_values('pval_adj')

    print(f'Kruskal-Wallis: {len(results_df)} interactions tested, '
          f'{(results_df["pval_adj"] < 0.05).sum() if len(results_df) else 0} significant at FDR < 0.05')
    return results_df


def test_one_vs_rest(interaction_matrix, clusters):
    """
    Two-sided Mann-Whitney U test per interaction, comparing each cluster
    against all remaining samples, with Benjamini-Hochberg correction within
    each cluster.
    """
    unique_clusters = sorted(clusters.unique())
    rows = []

    for cluster in unique_clusters:
        in_cluster = interaction_matrix.loc[clusters == cluster]
        out_cluster = interaction_matrix.loc[clusters != cluster]

        for feature in interaction_matrix.columns:
            a = in_cluster[feature].values
            b = out_cluster[feature].values
            try:
                stat, p = mannwhitneyu(a, b)
            except Exception:
                continue
            rows.append({
                'cluster': cluster,
                'feature': feature,
                'pval': p,
                'mean_in': a.mean(),
                'mean_out': b.mean(),
                'diff': a.mean() - b.mean()
            })

    ovr_df = pd.DataFrame(rows)
    if len(ovr_df):
        ovr_df['pval_adj'] = ovr_df.groupby('cluster')['pval'].transform(
            lambda x: multipletests(x, method='fdr_bh')[1]
        )
        ovr_df = ovr_df.sort_values(['cluster', 'diff'], ascending=[True, False])

    print(f'One-versus-rest: {len(ovr_df)} tests, '
          f'{(ovr_df["pval_adj"] < 0.05).sum() if len(ovr_df) else 0} significant at FDR < 0.05')
    return ovr_df


def simplify_celltype(name, rules):
    """
    Map detailed cell state names to shorter labels for display. Order matters:
    the first matching rule wins.
    """
    for pattern, new_name in rules:
        if re.search(pattern, name, flags=re.IGNORECASE):
            return new_name
    return name


def top_interactions_per_cluster(interaction_matrix, clusters, ovr_df, n_top, rename_rules=None):
    """
    Mean score per cluster for the n_top interactions with the largest
    one-versus-rest difference in each cluster.
    """
    top_features = (ovr_df
                    .sort_values(['cluster', 'diff'], ascending=[True, False])
                    .groupby('cluster')
                    .head(n_top)['feature'].unique())

    heatmap_data = (interaction_matrix[top_features]
                    .groupby(clusters)
                    .mean()
                    .T)

    if rename_rules:
        idx = heatmap_data.index
        sender = idx.str.split('|').str[0].str.split('__').str[0]
        receiver = idx.str.split('|').str[0].str.split('__').str[1]
        ligand = idx.str.split('|').str[1].str.split('__').str[0]
        receptor = idx.str.split('|').str[1].str.split('__').str[1]

        # simplify cell type names only, never ligand or receptor gene names
        sender = sender.map(lambda x: simplify_celltype(x, rename_rules))
        receiver = receiver.map(lambda x: simplify_celltype(x, rename_rules))

        heatmap_data.index = sender + '__' + receiver + '|' + ligand + '__' + receptor

    return heatmap_data


def plot_heatmap(heatmap_data, output_path):
    """
    Z-scored heatmap of the top discriminating interactions per cluster.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    heatmap_z = (heatmap_data
                 .subtract(heatmap_data.mean(axis=1), axis=0)
                 .divide(heatmap_data.std(axis=1).replace(0, 1), axis=0))

    fig, ax = plt.subplots(figsize=(10, max(6, len(heatmap_z) * 0.25)))
    sns.heatmap(heatmap_z, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                linewidths=0.3, cbar_kws={'label': 'Z-score'}, ax=ax)
    ax.set_title('Top discriminating interactions per cluster')
    ax.set_xlabel('Interaction group')
    ax.set_ylabel('')
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Group samples by their cell-cell communication landscape.')
    parser.add_argument('--adata', required=True,
                        help='Path to the .h5ad file')
    parser.add_argument('--output-dir', required=True,
                        help='Directory for output tables and figures')
    parser.add_argument('--liana-dir', default=None,
                        help='Directory of per-sample pickled LIANA results, used when '
                             'results are not stored in adata.uns')
    parser.add_argument('--sample-col', default='sample',
                        help='Column in adata.obs holding the sample identifier')
    parser.add_argument('--group', default='cell_states_detailed',
                        help='Cell type annotation level the LIANA results were computed at; '
                             'used to build the pickle filename')
    parser.add_argument('--uns-prefix', default=None,
                        help='Prefix of the adata.uns keys holding LIANA results. '
                             'Default: all_liana_<group>')
    parser.add_argument('--threshold', type=float, default=0.1,
                        help='Maximum magnitude and specificity rank for an interaction to be kept')
    parser.add_argument('--n-components', type=int, default=50,
                        help='Number of truncated SVD components to compute')
    parser.add_argument('--n-pcs', type=int, default=30,
                        help='Number of components used for the graph and the embedding')
    parser.add_argument('--n-neighbors', type=int, default=10,
                        help='Number of neighbours for the graph and for UMAP')
    parser.add_argument('--resolution', type=float, default=1.0,
                        help='Leiden resolution')
    parser.add_argument('--n-top', type=int, default=10,
                        help='Number of top interactions per cluster shown in the heatmap')
    parser.add_argument('--no-umap', action='store_true',
                        help='Skip the UMAP embedding')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    uns_prefix = args.uns_prefix if args.uns_prefix is not None else f'all_liana_{args.group}'

    print('-> Reading dataset')
    adata = sc.read_h5ad(args.adata)
    print(f'   {adata.shape[0]} cells x {adata.shape[1]} genes')

    print('-> Loading per-sample LIANA results')
    liana_results = load_liana_results(adata, args.sample_col, uns_prefix,
                                       args.liana_dir, args.group)

    print('-> Building interaction matrix')
    interaction_matrix = build_interaction_matrix(liana_results, args.threshold)
    interaction_matrix.to_csv(os.path.join(args.output_dir, 'interaction_matrix.csv'))

    print('-> Reducing dimensionality')
    components, svd = reduce_dimensions(interaction_matrix, args.n_components)

    print('-> Clustering samples')
    clusters = cluster_samples(components, args.n_pcs, args.n_neighbors, args.resolution)

    embedding = components.iloc[:, :args.n_pcs].copy()
    if not args.no_umap:
        print('-> Computing UMAP')
        umap_df = compute_umap(components, args.n_pcs, args.n_neighbors)
        embedding = embedding.join(umap_df)
    embedding['cluster'] = clusters
    embedding.to_csv(os.path.join(args.output_dir, 'sample_embedding.csv'))

    print('-> Testing interactions across clusters')
    kruskal_df = test_across_clusters(interaction_matrix, clusters)
    kruskal_df.to_csv(os.path.join(args.output_dir, 'kruskal_results.csv'), index=False)

    print('-> Testing interactions one versus rest')
    ovr_df = test_one_vs_rest(interaction_matrix, clusters)
    ovr_df.to_csv(os.path.join(args.output_dir, 'one_vs_rest_results.csv'), index=False)

    if len(ovr_df):
        print('-> Summarizing top interactions')
        heatmap_data = top_interactions_per_cluster(interaction_matrix, clusters,
                                                    ovr_df, args.n_top)
        heatmap_data.to_csv(os.path.join(args.output_dir, 'top_interactions_heatmap.csv'))
        plot_heatmap(heatmap_data, os.path.join(args.output_dir, 'top_interactions_heatmap.png'))

    print(f'Done. Results saved in {args.output_dir}')


if __name__ == '__main__':
    main()
