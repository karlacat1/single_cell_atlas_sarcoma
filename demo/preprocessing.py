"""
Quality control, filtering and normalization for a single sample.

Implements the per-sample preprocessing steps: initial gene and cell filtering,
QC metric calculation and plotting, outlier detection by median absolute deviation,
doublet detection with Scrublet, normalization and log transformation, highly
variable gene selection, cell cycle scoring, regression of technical covariates,
scaling and PCA.

External resource files are read from the directory given by the SARCOMA_RESOURCES
environment variable (default: ../resources relative to this file). See
resources/README.md for how to obtain them.
"""

import os

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib.pyplot import rc_context
from scipy.stats import median_abs_deviation as mad

RESOURCES_DIR = os.environ.get(
    'SARCOMA_RESOURCES',
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'resources')
)

RIBO_GENES_FILE = os.path.join(RESOURCES_DIR, 'KEGG_RIBOSOME.v2023.1.Hs.txt')
CELL_CYCLE_GENES_FILE = os.path.join(RESOURCES_DIR, 'regev_lab_cell_cycle_genes.txt')


class PrepProcessing:

    def __init__(self, adata, output_directory, filename):
        """
        :param adata: AnnData object for one sample
        :param output_directory: directory for QC plots
        :param filename: sample name, used as a prefix for output files
        """
        self.adata_prep = adata
        self.output_directory = output_directory
        self.filename = filename

        os.makedirs(output_directory, exist_ok=True)

    def initial_filtering(self, min_genes, min_cells):
        """
        Remove cells with fewer than min_genes detected genes and genes detected
        in fewer than min_cells cells.
        """
        print('---> Initial filtering')
        sc.pp.filter_cells(self.adata_prep, min_genes=min_genes)
        sc.pp.filter_genes(self.adata_prep, min_cells=min_cells)

    def highest_expr_genes(self, nr_genes, show):
        """
        Plot the genes with the highest fraction of counts.
        """
        print('---> Computing Top 20 highest expressed genes')
        with rc_context({'figure.figsize': (10, 8)}):
            sc.pl.highest_expr_genes(self.adata_prep, n_top=nr_genes, show=show)
        plt.title('Top 20 highly expressed genes')
        plt.savefig(self.output_directory + '/' + self.filename + '_top_20_highest_expr.png',
                    dpi=100, bbox_inches='tight')
        plt.close()

    def calculate_metrics(self):
        """
        Flag mitochondrial, ribosomal, hemoglobin and MALAT1 genes and compute
        per-cell QC metrics.
        """
        print('---> Calculating QC metrics ... ')
        ribo_genes = pd.read_table(RIBO_GENES_FILE, skiprows=2, header=None)
        self.adata_prep.var['malat'] = self.adata_prep.var_names.str.startswith('MALAT')
        self.adata_prep.var['mt'] = self.adata_prep.var_names.str.startswith('MT-')
        self.adata_prep.var['hb'] = self.adata_prep.var_names.str.contains("^HB[^(P)]")
        self.adata_prep.var['ribo'] = self.adata_prep.var_names.isin(ribo_genes[0].values)

        sc.pp.calculate_qc_metrics(self.adata_prep, qc_vars=["mt", "ribo", "hb"],
                                   inplace=True, percent_top=[20], log1p=True)

        remove = ['total_counts_mt', 'log1p_total_counts_mt', 'total_counts_ribo',
                  'log1p_total_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb']
        self.adata_prep.obs = self.adata_prep.obs[[x for x in self.adata_prep.obs.columns if x not in remove]]

    def qc_violin_plots(self, show=False, name_ending=''):
        """
        Violin plots of the main QC metrics.
        """
        print('---> Computing QC violin plots')
        with rc_context({'figure.figsize': (10, 8)}):
            sc.pl.violin(self.adata_prep,
                         ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
                         jitter=0.4, multi_panel=True, show=show)
        plt.savefig(self.output_directory + '/' + self.filename + '_violin_qc_plots_' + name_ending + '.png',
                    dpi=120, bbox_inches='tight')
        plt.close()

    def genes_vs_count_depth_plot(self):
        """
        Number of detected genes against count depth, colored by mitochondrial
        percentage. High mitochondrial fractions are expected mainly in low-count
        cells with few detected genes.
        """
        print('---> Plotting Gene VS count depth')
        with rc_context({'figure.figsize': (7, 7)}):
            sc.pl.scatter(self.adata_prep, x='total_counts', y='n_genes_by_counts',
                          title='Nr Genes vs Count depth',
                          color='pct_counts_mt', show=False)
        plt.savefig(self.output_directory + '/' + self.filename + '_genes_vs_count_depth.png',
                    bbox_inches='tight')

        print('Correlation between count depth and nr of genes per cell: ',
              stats.pearsonr(self.adata_prep.obs['total_counts'],
                             self.adata_prep.obs['n_genes_by_counts'])[0])

    def mt_vs_count_depth_plot(self, mt_cutoff):
        """
        Mitochondrial percentage against count depth, with the minimum
        mitochondrial threshold drawn for reference.
        """
        print('---> Plotting MT percentage VS count depth')
        with rc_context({'figure.figsize': (6, 6)}):
            plt.scatter(self.adata_prep.obs['total_counts'], self.adata_prep.obs['pct_counts_mt'],
                        color='grey', s=2)
            if mt_cutoff > 0:
                plt.axhline(mt_cutoff, color='r', ls='--', label='possible cutoff')
                plt.legend()
            plt.title('Percentage of MT vs Count Depth')
            plt.ylabel('pct_counts_mt')
            plt.xlabel('total_counts')

        plt.savefig(self.output_directory + '/' + self.filename + '_mt_vs_count_depth.png',
                    bbox_inches='tight')
        print('Pearson Correlation count depth vs percentage of MT: ',
              stats.pearsonr(self.adata_prep.obs['total_counts'],
                             self.adata_prep.obs['pct_counts_mt'])[0])

    def outlier_mads(self, mt_cutoff):
        """
        Flag cells deviating by more than 5 median absolute deviations from the
        sample median in log1p total counts, log1p number of detected genes,
        percentage of counts in the top 20 genes, or mitochondrial percentage.

        Mitochondrial percentage is treated as upper-tail only. Because many
        nuclei have almost no mitochondrial counts, the MAD-derived threshold can
        fall very low; mt_cutoff is used as a floor so that cells below this
        percentage are never flagged on this criterion alone.
        """

        def mad_outlier(adata, metric, nmads, upper_only=False):

            M = adata.obs[metric]

            if not upper_only:
                return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

            threshold = np.median(M) + nmads * mad(M)
            if threshold < mt_cutoff:
                print('VERY LOW MADS threshold:', threshold)
                threshold = mt_cutoff

            return (M > threshold)

        self.adata_prep.obs["mad_outlier"] = (
                mad_outlier(self.adata_prep, "log1p_total_counts", 5)
                | mad_outlier(self.adata_prep, "log1p_n_genes_by_counts", 5)
                | mad_outlier(self.adata_prep, "pct_counts_in_top_20_genes", 5)
                | mad_outlier(self.adata_prep, "pct_counts_mt", 5, upper_only=True)
        )

    def doublet_detection_and_removal(self, show=False):
        """
        Score doublets with Scrublet (expected doublet rate 0.07). Cells are
        flagged here and removed in filtering_cells_qc_metrics.
        """
        print('---> Doublet detection:')
        scrub = scr.Scrublet(self.adata_prep.X, expected_doublet_rate=0.07)
        out = scrub.scrub_doublets()

        dfg = pd.DataFrame({'doublet_score': out[0], 'predicted_doublets': out[1]},
                           index=self.adata_prep.obs.index)
        print(dfg.predicted_doublets.sum(), " predicted doublets")
        self.adata_prep.obs['doublet_scores'] = dfg['doublet_score']
        self.adata_prep.obs['predicted_doublets'] = dfg['predicted_doublets'].astype(str)

        with rc_context({'figure.figsize': (10, 8)}):
            sc.pl.violin(self.adata_prep, 'n_genes_by_counts',
                         jitter=0.4, groupby='predicted_doublets', rotation=45,
                         show=show,)
        plt.savefig(self.output_directory + '/' + self.filename + '_doublet_compare.png',
                    dpi=120, bbox_inches='tight')
        plt.close()

        scrub.plot_histogram()
        plt.savefig(self.output_directory + '/' + self.filename + '_doublet_histogram.png',
                    dpi=120, bbox_inches='tight')
        plt.close()

    def filtering_cells_qc_metrics(self):
        """
        Apply the QC filters: remove cells above 8% mitochondrial counts, cells
        flagged as MAD outliers, and predicted doublets; drop mitochondrial and
        MALAT1 genes from the feature space.
        """
        print('---> Filtering part 2:')
        self.adata_prep = self.adata_prep[self.adata_prep.obs.pct_counts_mt < 8]

        self.adata_prep = self.adata_prep[~self.adata_prep.obs['mad_outlier'], :]

        self.adata_prep = self.adata_prep[:, -self.adata_prep.var['mt']]
        self.adata_prep = self.adata_prep[:, -self.adata_prep.var['malat']]

        self.adata_prep = self.adata_prep[self.adata_prep.obs['predicted_doublets'] == 'False', :]

    def norm_and_log(self):
        """
        Normalize each cell to 10,000 counts and log1p-transform.
        """
        print('---> Normalizing and log-transforming the data ...')
        sc.pp.normalize_total(self.adata_prep, target_sum=1e4)
        sc.pp.log1p(self.adata_prep)
        return self.adata_prep

    def compute_highly_variable(self):
        """
        Identify highly variable genes (Scanpy default, flavor='seurat').
        """
        sc.pp.highly_variable_genes(self.adata_prep)

    def cell_cycle_score(self):
        """
        Score S and G2/M phase using the gene sets of Tirosh et al. (2016).
        The first 43 genes in the file are S phase, the remainder G2/M.
        """
        print('---> Computing Cell cycle score ')
        cell_cycle_genes = [x.strip() for x in open(CELL_CYCLE_GENES_FILE)]
        s_genes = cell_cycle_genes[:43]
        g2m_genes = cell_cycle_genes[43:]
        scaled_data = sc.pp.scale(self.adata_prep, copy=True)

        sc.tl.score_genes_cell_cycle(scaled_data, s_genes=s_genes, g2m_genes=g2m_genes)

        self.adata_prep.obs['S_score'] = scaled_data.obs['S_score']
        self.adata_prep.obs['G2M_score'] = scaled_data.obs['G2M_score']
        self.adata_prep.obs['phase'] = scaled_data.obs['phase']

    def filter_hvariable_regress(self):
        """
        Subset to highly variable genes and regress out total counts and
        mitochondrial percentage.
        """
        print('---> Filtering highly variable ')
        self.adata_prep = self.adata_prep[:, self.adata_prep.var.highly_variable]
        print('---> Regress out effect of total counts per cell.. ')
        sc.pp.regress_out(self.adata_prep, ['total_counts', 'pct_counts_mt'])
        return self.adata_prep

    def scale(self, max_deviation):
        """
        Scale to unit variance, clipping at max_deviation standard deviations.
        """
        print('---> Scaling Data')
        sc.pp.scale(self.adata_prep, max_value=max_deviation)

    def pca_inspection(self):
        """
        Principal component analysis on the scaled data.
        """
        print('---> Applying PCA')
        sc.tl.pca(self.adata_prep, svd_solver='arpack', random_state=42)
