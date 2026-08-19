"""
Group samples by their cell-cell communication landscape.

Usage:
  # per sample
  python interaction_grouping.py --adata atlas.h5ad --liana-dir ccc/per_sample \
      --group sample --split-by sample --output-dir results/

  # per sarcoma entity
  python interaction_grouping.py --adata atlas.h5ad --liana-dir ccc/per_entity \
      --group cell_types --split-by entity --sample-col entity --output-dir results/

  # whole cohort
  python interaction_grouping.py --adata atlas.h5ad --liana-dir ccc/all \
      --group cell_types --output-dir results/

Note:
  Whole-cohort LIANA produces one interaction table for the entire cohort.
  The clustering workflow below requires multiple independent observations
  (e.g. samples/entities), so whole-cohort mode is loaded and reported but
  cannot perform sample-level clustering.
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


def load_liana_results(adata, sample_col, liana_dir, group, split_by=None):
    """Load whole-cohort or split LIANA results according to --split-by."""

    if liana_dir is None:
        raise ValueError('--liana-dir is required when loading LIANA pickle results.')

    # Whole cohort: run_liana.py --split-by omitted
    if split_by is None:
        path = os.path.join(liana_dir, f'uns_liana_{group}.pkl')

        if not os.path.exists(path):
            raise FileNotFoundError(f'Whole-cohort LIANA result not found: {path}')

        with open(path, 'rb') as f:
            df = pd.DataFrame(pickle.load(f)).copy()

        missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(f'Whole-cohort LIANA results are missing columns: {missing}')

        print(f'Loaded whole-cohort LIANA results: {len(df)} interactions')
        return {'all': df}

    # Split mode: run_liana.py --split-by <column>
    if split_by not in adata.obs.columns:
        raise ValueError(f'--split-by column not found in adata.obs: {split_by}')

    levels = adata.obs[split_by].astype(str).unique()
    results = {}

    for level in levels:
        path = os.path.join(liana_dir, str(level), f'uns_liana_{level}_{group}.pkl')

        if os.path.exists(path):
            with open(path, 'rb') as f:
                results[str(level)] = pd.DataFrame(pickle.load(f)).copy()
        else:
            print(f'No LIANA results for {level}')

    if not results:
        raise ValueError(f'No LIANA results found in {liana_dir}')

    for level, df in results.items():
        missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(f'LIANA results for {level} are missing columns: {missing}')

    print(f'Loaded LIANA results for {len(results)} of {len(levels)} {split_by} levels')
    return results


def build_interaction_matrix(liana_results, threshold):
    """Build the observation-by-interaction score matrix."""

    all_rows = []

    for sample, liana_res in liana_results.items():
        liana_res = liana_res.copy()

        liana_res['feature'] = (liana_res['source'].astype(str) + '__' +
                                liana_res['target'].astype(str) + '|' +
                                liana_res['ligand_complex'].astype(str) + '__' +
                                liana_res['receptor_complex'].astype(str))

        keep = ((liana_res['magnitude_rank'] <= threshold) &
                (liana_res['specificity_rank'] <= threshold))
        liana_res = liana_res[keep]

        liana_res['score'] = ((1 - liana_res['magnitude_rank']) *
                              (1 - liana_res['specificity_rank']))
        liana_res['sample'] = sample

        all_rows.append(liana_res[['sample', 'feature', 'score']])

    if not all_rows:
        raise ValueError(f'No interactions passed threshold={threshold}.')

    long_df = pd.concat(all_rows, ignore_index=True)

    interaction_matrix = long_df.pivot_table(
        index='sample', columns='feature', values='score',
        aggfunc='max', fill_value=0
    )

    if interaction_matrix.shape[0] > 1:
        nonzero_var = interaction_matrix.columns[interaction_matrix.var(axis=0) > 0]
        interaction_matrix = interaction_matrix[nonzero_var]

    print(f'Interaction matrix: {interaction_matrix.shape[0]} observations x '
          f'{interaction_matrix.shape[1]} interactions')

    return interaction_matrix


def reduce_dimensions(interaction_matrix, n_components, random_state=42):
    """Scale without mean-centering and reduce with truncated SVD."""

    if interaction_matrix.shape[0] < 2:
        raise ValueError('At least 2 observations are required for dimensionality reduction. '
                         'Whole-cohort LIANA produces only one observation.')

    if interaction_matrix.shape[1] < 1:
        raise ValueError('No variable interactions remain after filtering. '
                         'Try increasing --threshold.')

    scaler = StandardScaler(with_mean=False)
    x_scaled = scaler.fit_transform(interaction_matrix.values)

    max_components = min(interaction_matrix.shape) - 1

    if max_components < 1:
        raise ValueError('Not enough observations/features for SVD.')

    n_components = min(n_components, max_components)

    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    coords = svd.fit_transform(x_scaled)

    components = pd.DataFrame(
        coords, index=interaction_matrix.index,
        columns=[f'PC{i + 1}' for i in range(n_components)]
    )

    print(f'Truncated SVD: {n_components} components, '
          f'{svd.explained_variance_ratio_.sum():.1%} of variance explained')

    return components, svd


def cluster_samples(components, n_pcs, n_neighbors, resolution, random_state=42):
    """Cluster observations using Leiden."""

    if components.shape[0] < 3:
        raise ValueError('At least 3 observations are recommended for clustering.')

    n_pcs = min(n_pcs, components.shape[1])
    n_neighbors = min(n_neighbors, components.shape[0] - 1)

    tmp = ad.AnnData(X=components.iloc[:, :n_pcs].values)
    tmp.obs.index = components.index

    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, use_rep='X',
                    random_state=random_state)
    sc.tl.leiden(tmp, resolution=resolution, random_state=random_state)

    clusters = tmp.obs['leiden'].rename('cluster')
    clusters.index = components.index

    print(f'Leiden clustering: {clusters.nunique()} interaction groups')
    print(clusters.value_counts().sort_index().to_string())

    return clusters


def compute_umap(components, n_pcs, n_neighbors, min_dist=0.1, random_state=42):
    """Compute UMAP embedding."""

    import umap

    n_pcs = min(n_pcs, components.shape[1])
    n_neighbors = min(n_neighbors, components.shape[0] - 1)

    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist,
                        n_components=2, random_state=random_state)

    coords = reducer.fit_transform(components.iloc[:, :n_pcs].values)

    return pd.DataFrame(coords, index=components.index,
                        columns=['UMAP1', 'UMAP2'])


def test_across_clusters(interaction_matrix, clusters):
    """Kruskal-Wallis test per interaction across clusters."""

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

    results_df = pd.DataFrame(results)

    if len(results_df):
        results_df['pval_adj'] = multipletests(
            results_df['pval'], method='fdr_bh'
        )[1]
        results_df = results_df.sort_values('pval_adj')

    print(f'Kruskal-Wallis: {len(results_df)} interactions tested, '
          f'{(results_df["pval_adj"] < 0.05).sum() if len(results_df) else 0} '
          f'significant at FDR < 0.05')

    return results_df


def test_one_vs_rest(interaction_matrix, clusters):
    """Mann-Whitney U test for each cluster versus the rest."""

    unique_clusters = sorted(clusters.unique())
    rows = []

    for cluster in unique_clusters:
        in_cluster = interaction_matrix.loc[clusters == cluster]
        out_cluster = interaction_matrix.loc[clusters != cluster]

        for feature in interaction_matrix.columns:
            a = in_cluster[feature].values
            b = out_cluster[feature].values

            if len(a) < 2 or len(b) < 2:
                continue

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
        ovr_df = ovr_df.sort_values(
            ['cluster', 'diff'], ascending=[True, False]
        )

    print(f'One-versus-rest: {len(ovr_df)} tests, '
          f'{(ovr_df["pval_adj"] < 0.05).sum() if len(ovr_df) else 0} '
          f'significant at FDR < 0.05')

    return ovr_df


def simplify_celltype(name, rules):
    """Map detailed cell state names to shorter labels."""

    for pattern, new_name in rules:
        if re.search(pattern, name, flags=re.IGNORECASE):
            return new_name
    return name


def top_interactions_per_cluster(interaction_matrix, clusters, ovr_df,
                                  n_top, rename_rules=None):
    """Summarize top interactions per cluster."""

    top_features = (ovr_df.sort_values(
        ['cluster', 'diff'], ascending=[True, False]
    ).groupby('cluster').head(n_top)['feature'].unique())

    heatmap_data = interaction_matrix[top_features].groupby(clusters).mean().T

    if rename_rules:
        idx = heatmap_data.index
        sender = idx.str.split('|').str[0].str.split('__').str[0]
        receiver = idx.str.split('|').str[0].str.split('__').str[1]
        ligand = idx.str.split('|').str[1].str.split('__').str[0]
        receptor = idx.str.split('|').str[1].str.split('__').str[1]

        sender = sender.map(lambda x: simplify_celltype(x, rename_rules))
        receiver = receiver.map(lambda x: simplify_celltype(x, rename_rules))

        heatmap_data.index = (
            sender + '__' + receiver + '|' + ligand + '__' + receptor
        )

    return heatmap_data


def plot_heatmap(heatmap_data, output_path):
    """Plot z-scored heatmap."""

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    heatmap_z = (heatmap_data.subtract(
        heatmap_data.mean(axis=1), axis=0
    ).divide(
        heatmap_data.std(axis=1).replace(0, 1), axis=0
    ))

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

    parser.add_argument('--adata', required=True, help='Path to the .h5ad file')
    parser.add_argument('--output-dir', required=True, help='Directory for output tables and figures')
    parser.add_argument('--liana-dir', required=True, help='Directory containing LIANA pickle results')
    parser.add_argument('--sample-col', default='sample', help='Column in adata.obs holding the sample identifier')
    parser.add_argument('--group', default='cell_states_detailed', help='Annotation level used for LIANA')
    parser.add_argument('--threshold', type=float, default=0.1, help='Maximum magnitude and specificity rank')
    parser.add_argument('--n-components', type=int, default=50, help='Number of truncated SVD components')
    parser.add_argument('--n-pcs', type=int, default=30, help='Number of components used for clustering/UMAP')
    parser.add_argument('--n-neighbors', type=int, default=10, help='Number of neighbours')
    parser.add_argument('--resolution', type=float, default=1.0, help='Leiden resolution')
    parser.add_argument('--n-top', type=int, default=10, help='Number of top interactions per cluster')
    parser.add_argument('--no-umap', action='store_true', help='Skip UMAP')
    parser.add_argument('--split-by', default=None, help='Column used to split LIANA analysis; omit for whole cohort')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print('-> Reading dataset')
    adata = sc.read_h5ad(args.adata)
    print(f'   {adata.shape[0]} cells x {adata.shape[1]} genes')

    if args.split_by is None:
        print('-> Loading whole-cohort LIANA results')
    else:
        print(f'-> Loading LIANA results split by {args.split_by}')

    liana_results = load_liana_results(
        adata, args.sample_col, args.liana_dir, args.group, args.split_by
    )

    print('-> Building interaction matrix')
    interaction_matrix = build_interaction_matrix(liana_results, args.threshold)
    interaction_matrix.to_csv(os.path.join(args.output_dir, 'interaction_matrix.csv'))

    # Whole-cohort LIANA has one observation, so there is nothing to cluster.
    if args.split_by is None:
        print('-> Whole-cohort LIANA loaded successfully.')
        print('   Skipping clustering/statistical analysis because the whole-cohort '
              'result contains one observation.')
        print(f'Done. Result saved in {args.output_dir}')
        return

    print('-> Reducing dimensionality')
    components, svd = reduce_dimensions(interaction_matrix, args.n_components)

    print('-> Clustering samples')
    clusters = cluster_samples(
        components, args.n_pcs, args.n_neighbors, args.resolution
    )

    embedding = components.iloc[:, :min(args.n_pcs, components.shape[1])].copy()

    if not args.no_umap:
        print('-> Computing UMAP')
        umap_df = compute_umap(
            components, args.n_pcs, args.n_neighbors
        )
        embedding = embedding.join(umap_df)

    embedding['cluster'] = clusters
    embedding.to_csv(
        os.path.join(args.output_dir, 'sample_embedding.csv')
    )

    print('-> Testing interactions across clusters')
    kruskal_df = test_across_clusters(interaction_matrix, clusters)
    kruskal_df.to_csv(
        os.path.join(args.output_dir, 'kruskal_results.csv'),
        index=False
    )

    print('-> Testing interactions one versus rest')
    ovr_df = test_one_vs_rest(interaction_matrix, clusters)
    ovr_df.to_csv(
        os.path.join(args.output_dir, 'one_vs_rest_results.csv'),
        index=False
    )

    if len(ovr_df):
        print('-> Summarizing top interactions')
        heatmap_data = top_interactions_per_cluster(
            interaction_matrix, clusters, ovr_df, args.n_top
        )
        heatmap_data.to_csv(
            os.path.join(args.output_dir, 'top_interactions_heatmap.csv')
        )
        plot_heatmap(
            heatmap_data,
            os.path.join(args.output_dir, 'top_interactions_heatmap.png')
        )

    print(f'Done. Results saved in {args.output_dir}')


if __name__ == '__main__':
    main()
