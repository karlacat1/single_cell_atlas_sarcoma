"""
Dimensionality reduction and plotting for a single sample.

Wraps the PCA, neighborhood graph, UMAP and QC visualization steps used by the
per-sample pipeline. All figures are written to the output directory with the
sample name as a prefix.
"""

import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context


class Visualization:
    def __init__(self, adata, output_directory, filename):
        """
        :param adata: preprocessed AnnData object
        :param output_directory: directory for figures
        :param filename: sample name, used as a prefix for output files
        """
        sc.settings.set_figure_params(dpi=200, facecolor='white')
        self.adata_vis = adata
        self.output_directory = output_directory
        self.filename = filename
        # Scanpy stores the log base in .uns; reset it so that downstream
        # differential expression does not rescale log fold changes.
        self.adata_vis.uns['log1p']['base'] = None

    def pca_inspection_plot(self):
        """
        Elbow plot of the PCA variance ratio, used to choose the number of
        principal components for the neighborhood graph.
        """
        print('---> Plotting PCA')
        sc.pl.pca_variance_ratio(self.adata_vis, log=True, show=False, n_pcs=60)
        plt.savefig(self.output_directory + '/' + self.filename + '_pca_elbow.png', bbox_inches='tight')

    def pca_plot(self, clustering, method):
        """
        Scatter plot of the first two principal components.

        :param clustering: whether to color by a clustering result
        :param method: name of the obs column holding the clustering
        """
        print('--->Plotting PCA plot')
        if clustering and method != "":
            sc.pl.pca(self.adata_vis, show=False, color=method, title='PCA plot for ' + self.filename)
            plt.savefig(self.output_directory + '/' + self.filename + '_pca_clustering.png', bbox_inches='tight')
        else:
            sc.pl.pca(self.adata_vis, show=False, title='PCA plot for ' + self.filename)
            plt.savefig(self.output_directory + '/' + self.filename + '_pca.png', bbox_inches='tight')

    def knn_graph(self, npcs, n_neighbors):
        """
        Build the nearest-neighbor graph on the first npcs principal components.
        """
        print('--->Computing neighbour graph')
        sc.pp.neighbors(self.adata_vis, n_neighbors=n_neighbors, use_rep='X_pca', n_pcs=npcs, random_state=42)

    def umap_alg(self):
        """
        Compute the UMAP embedding from the neighborhood graph.
        """
        sc.tl.umap(self.adata_vis, random_state=42)

    def umap_visualization(self, color, title, palette=None):
        """
        Plot the UMAP embedding colored by one or several obs columns or genes.

        :param color: obs column or gene name, or a list of them
        :param title: plot title, or a list of titles matching `color`
        :param palette: optional color palette
        """
        print('---> Visualizing data using umap...')
        out_file = self.output_directory + '/' + self.filename

        if type(title) != str:
            sc.pl.umap(self.adata_vis, color=color, title=title, frameon=False, show=False,
                       legend_fontsize=10, legend_fontoutline=2, palette=palette)
            plt.savefig(out_file + '_multiple_umap.png', bbox_inches='tight')
        elif palette is None:
            sc.pl.umap(self.adata_vis, color=color, title=title, frameon=False, show=False,
                       legend_fontsize=10, legend_fontoutline=2)
            plt.savefig(out_file + '_' + title + '_umap2.png', bbox_inches='tight')
        else:
            sc.pl.umap(self.adata_vis, color=color, title=title, frameon=False, show=False,
                       legend_fontsize=10, legend_fontoutline=2, palette=palette)
            plt.savefig(out_file + '_' + title + '_umap.png', bbox_inches='tight')

    def cell_cycle_inspection(self, clustering_done):
        """
        Violin plots of the S and G2/M phase scores, grouped by cluster if
        clustering has been performed.
        """
        with rc_context({'figure.figsize': (10, 8)}):
            if clustering_done:
                sc.pl.violin(self.adata_vis, ['S_score', 'G2M_score'],
                             jitter=0.4, groupby='leiden', rotation=45, show=False)
            else:
                sc.pl.violin(self.adata_vis, ['S_score', 'G2M_score'],
                             jitter=0.4, rotation=45, show=False)
        plt.savefig(self.output_directory + '/' + self.filename + '_cell_cycle_score.png',
                    dpi=100, bbox_inches='tight')
        plt.close()
