"""
Leiden clustering of a preprocessed sample.
"""

import numpy as np
import scanpy as sc


class Clustering:
    def __init__(self, adata, output_directory, filename):
        """
        :param adata: preprocessed AnnData object with a neighborhood graph
        :param output_directory: directory for output files
        :param filename: sample name, used as a prefix for output files
        """
        self.expression_data = adata
        self.output_directory = output_directory
        self.filename = filename

    def clustering_leiden(self, res):
        """
        Cluster cells with the Leiden algorithm at the given resolution.
        """
        print('---> Applying leiden clustering...')
        sc.tl.leiden(self.expression_data, random_state=42, resolution=res)
        print('Nr of clusters found: {}'.format(len(np.unique(self.expression_data.obs['leiden']))))
