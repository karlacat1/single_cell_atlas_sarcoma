"""
Define TME archetypes by hierarchical clustering of non-malignant composition.

Each sample is described by the proportion of each non-malignant cell type among
all its non-malignant cells. Samples are clustered with Ward linkage on Euclidean
distances and cut into a fixed number of clusters.

Input: .h5ad file whose .obs holds the sample identifier, the cell type
annotation and any sample-level columns to display as annotation tracks.

Output: <prefix>_archetypes.csv (assignment per sample),
        <prefix>_proportions.csv and <prefix>_archetypes.pdf in --output-dir.

Usage:
  python archetype_clustering.py --adata atlas.h5ad --output-dir results/ \
      --n-clusters 12 --annotation-cols Entity Stage_small Location Spatial
"""

import argparse
import itertools
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import to_hex
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list, linkage


def composition_matrix(obs, sample_col, annotation_col, malignant_label):
    """Proportion of each non-malignant cell type per sample."""
    obs = obs.astype(str).copy()
    obs['cells'] = obs.index.values
    non_malignant = obs[obs[annotation_col] != malignant_label]

    df = (non_malignant.groupby([sample_col, annotation_col]).count()['cells'] /
          non_malignant.groupby([sample_col])['cells'].count()).unstack().fillna(0)

    print(f'{df.shape[0]} samples x {df.shape[1]} non-malignant cell types')
    return df, non_malignant


def cluster_samples(df, n_clusters):
    """Ward linkage on Euclidean distances, cut into n_clusters."""
    linkage_matrix = linkage(df, method='ward', metric='euclidean')
    cluster_labels = fcluster(linkage_matrix, t=n_clusters, criterion='maxclust')
    print(f'Using fixed number of clusters: {n_clusters}')
    return linkage_matrix, cluster_labels


def make_color_map(values, cmap='tab20'):
    """One color per unique value, with empty strings shown in grey."""
    uniques = [v for v in pd.unique(values)]
    palette = plt.get_cmap(cmap).colors
    return {v: ('gainsboro' if v == '' else palette[i % len(palette)])
            for i, v in enumerate(uniques)}


def plot_archetypes(df, linkage_matrix, cluster_labels, annotations, cluster_colors,
                    celltype_colors, output_path, title):
    """Dendrogram, stacked composition bars and sample annotation tracks."""
    n = len(cluster_labels)

    # a dendrogram branch is colored when all its leaves share one cluster
    node_clusters = {i: {cluster_labels[i]} for i in range(n)}
    for i, (c1, c2, _, _) in enumerate(linkage_matrix):
        node_clusters[n + i] = node_clusters[int(c1)] | node_clusters[int(c2)]

    def link_color_func(node_id):
        clusters = node_clusters[node_id]
        return cluster_colors[list(clusters)[0] - 1] if len(clusters) == 1 else 'black'

    fig, ax = plt.subplots(1, 2, figsize=(20, 25), constrained_layout=True,
                           gridspec_kw={'width_ratios': [1, 4]})
    fig.subplots_adjust(wspace=0)

    dendrogram(linkage_matrix, orientation='left', ax=ax[0], labels=df.index,
               link_color_func=link_color_func)
    ax[0].axis('off')

    # reorder samples to match the dendrogram
    ax_bar = ax[1]
    cluster_order = leaves_list(linkage_matrix)
    df = df.iloc[cluster_order]
    cluster_labels_ordered = cluster_labels[cluster_order]

    df.plot(kind='barh', stacked=True, ax=ax_bar, color=celltype_colors,
            width=0.9, fontsize=14)
    ax_bar.grid(False, which='both', axis='both')
    ax_bar.legend(bbox_to_anchor=(1.1, 1.0), fontsize=14)
    ax_bar.set_ylabel('')

    # archetype track: one block per contiguous run of samples in the same cluster
    start_y = 0
    for label, group in itertools.groupby(cluster_labels_ordered):
        count = len(list(group))
        ax_bar.barh(left=1.01, y=[start_y + count / 2 - 0.5], width=0.05,
                    height=count, align='center', color=cluster_colors[label - 1])
        start_y += count

    # sample annotation tracks
    start_x = 1.08
    for column, (value_map, color_map) in annotations.items():
        for start_y, sample in enumerate(df.index):
            value = value_map.get(sample, '')
            ax_bar.barh(left=start_x, y=start_y, color=color_map.get(value, 'gainsboro'),
                        width=0.05, height=0.9, align='center')
        start_x += 0.05

    plt.title(title)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Define TME archetypes by hierarchical clustering.')
    parser.add_argument('--adata', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--prefix', default='archetypes')
    parser.add_argument('--n-clusters', type=int, default=12)
    parser.add_argument('--sample-col', default='sample')
    parser.add_argument('--annotation-col', default='cell_types',
                        help='Column holding the cell type annotation')
    parser.add_argument('--malignant-label', default='malignant',
                        help='Label of the malignant compartment, excluded from the composition')
    parser.add_argument('--annotation-cols', nargs='*', default=[],
                        help='Sample-level obs columns to display as annotation tracks')
    parser.add_argument('--celltype-colors', default=None,
                        help='Optional CSV with columns cell_type and color')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    df, non_malignant = composition_matrix(adata.obs, args.sample_col,
                                           args.annotation_col, args.malignant_label)
    df.to_csv(os.path.join(args.output_dir, f'{args.prefix}_proportions.csv'))

    linkage_matrix, cluster_labels = cluster_samples(df, args.n_clusters)

    assignments = pd.DataFrame({'sample': df.index, 'archetype': cluster_labels})
    assignments.to_csv(os.path.join(args.output_dir, f'{args.prefix}_archetypes.csv'), index=False)
    print(assignments['archetype'].value_counts().sort_index().to_string())

    # one color per archetype and per cell type
    cluster_colors = [to_hex(c) for c in
                      plt.cm.tab20(np.linspace(0, 1, len(np.unique(cluster_labels))))]

    celltype_colors = None
    if args.celltype_colors is not None:
        color_df = pd.read_csv(args.celltype_colors)
        celltype_colors = [dict(zip(color_df['cell_type'], color_df['color'])).get(c)
                           for c in df.columns]

    # sample-level annotation tracks
    annotations = {}
    for column in args.annotation_cols:
        if column not in non_malignant.columns:
            print(f'Annotation column not found, skipping: {column}')
            continue
        value_map = (non_malignant[[args.sample_col, column]].drop_duplicates()
                     .set_index(args.sample_col)[column].to_dict())
        annotations[column] = (value_map, make_color_map(list(value_map.values())))

    plot_archetypes(df, linkage_matrix, cluster_labels, annotations, cluster_colors,
                    celltype_colors,
                    os.path.join(args.output_dir, f'{args.prefix}_{args.n_clusters}_cluster.pdf'),
                    f'Clustering using fixed number of clusters: {args.n_clusters}')

    print(f'Done. Results saved in {args.output_dir}')


if __name__ == '__main__':
    main()
