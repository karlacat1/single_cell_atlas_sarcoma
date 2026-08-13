"""
Per-sample preprocessing, clustering and initial cell type annotation.

Runs the standard pipeline for a single sample: quality control and filtering,
doublet removal, normalization, feature selection, PCA, neighborhood graph,
Leiden clustering, UMAP, and an initial automated cell type annotation combining
differential expression, CellTypist, PanglaoDB and gene set enrichment into a
consensus call.

Input is either a Cell Ranger MEX directory or a preprocessed .h5ad file. The automated annotation produced here is a starting point for manual curation,
not the final annotation used in the manuscript. 

Usage:
    python main.py -d <path> [-n sample_name] [-o output_dir] [-j params.json] [-pr] [-dw]

Arguments:
    -d   Path to a Cell Ranger MEX directory (filtered_feature_*) or an .h5ad file.
    -j   Optional JSON file with parameter overrides. Defaults are given below;
         per-sample values (e.g. n_pcs, chosen from the PCA elbow plot) were
         supplied this way.
    -n   Sample name used as a prefix for output files. Defaults to the name of
         the input file or directory.
    -o   Output directory. Defaults to <sample name>_data_analysis_<date> created
         next to the input.
    -pr  Preprocessing only.
    -dw  Downstream analysis only; input must be a preprocessed .h5ad.
"""

import scanpy as sc
import warnings

warnings.simplefilter("ignore", category=FutureWarning)
import argparse
import os
from datetime import date
import time
import json
import preprocessing as pr
import visualization
import clustering
import cell_annotation
import utils

parameters = {
    "bins_count_depth": "150",
    "nr_bins_gene_depth": "150",
    "min_genes": "200",
    "min_cells": "3",
    "mt_cutoff": "5",
    "rib_cutoff": "8",
    "max_deviation": "10",
    "n_neighbours": "15",
    "n_pcs": "15",
    "leiden_resolution": "0.5"
}


def preprocessing(adata, output_directory, filename, nr_genes, nr_cells):
    """
    Quality control, filtering, normalization, feature selection and PCA.
    """

    print('-> Starting preprocessing ...')

    # 1) make variable names unique
    adata.var_names_make_unique()

    prep = pr.PrepProcessing(adata, output_directory + '/1.preprocessing', filename)

    # 2) initial filtering
    prep.initial_filtering(min_genes=int(parameters['min_genes']), min_cells=int(parameters['min_cells']))

    # 3) plot highest_expr_genes
    prep.highest_expr_genes(20, False)

    # 4) calculate QC metrics before filtering
    prep.calculate_metrics()

    # 5) plot QC metrics
    prep.qc_violin_plots(name_ending='initial')
    prep.mt_vs_count_depth_plot(int(parameters['mt_cutoff']))
    prep.genes_vs_count_depth_plot()

    # 6) flag outliers by median absolute deviation
    prep.outlier_mads(int(parameters['mt_cutoff']))

    # 7) doublet removal with Scrublet
    prep.doublet_detection_and_removal()
    prep.qc_violin_plots(name_ending='after_doublet')

    # 8) QC filtering
    prep.filtering_cells_qc_metrics()
    prep.qc_violin_plots(name_ending='after_filtering')

    # 9) store non-log-transformed counts for CNV inference
    utils.save_adata_h5ad(prep.adata_prep, output_directory, filename, 'non_log')
    utils.save_adata_for_r_script(adata, output_directory)

    # 10) normalization and log transformation
    adata = prep.norm_and_log()

    # 11) keep the log-normalized matrix as .raw for marker visualization
    adata.raw = adata

    # 12) identify highly variable genes
    prep.compute_highly_variable()

    # 13) cell cycle scoring
    prep.cell_cycle_score()

    # 14) subset to highly variable genes and regress out total counts and mitochondrial percentage
    adata = prep.filter_hvariable_regress()

    # 15) scale
    prep.scale(int(parameters['max_deviation']))

    print('Nr genes before filtering: ', nr_genes)
    print('Nr cells before filtering: ', nr_cells)
    nr_cells, nr_genes = adata.shape[0], adata.shape[1]
    print('Nr genes after filtering: ', nr_genes)
    print('Nr cells after filtering: ', nr_cells)

    # 16) PCA
    prep.pca_inspection()
    vis = visualization.Visualization(adata, output_directory + '/1.preprocessing', filename)
    vis.pca_inspection_plot()

    # 17) save preprocessed data
    utils.save_adata_h5ad(prep.adata_prep, output_directory, filename, 'preprocessed')

    return prep.adata_prep


def downstream_analysis(adata, output_directory, filename):
    """
    Neighborhood graph, Leiden clustering and UMAP embedding.
    """
    print('-> Starting initial DW analysis ...')
    os.makedirs(output_directory + '/2.downstream_analysis', exist_ok=True)

    # 1) neighborhood graph and Leiden clustering
    vis = visualization.Visualization(adata, output_directory + '/2.downstream_analysis', filename)
    vis.knn_graph(int(parameters['n_pcs']), int(parameters['n_neighbours']))
    cl = clustering.Clustering(adata, output_directory, filename)
    cl.clustering_leiden(float(parameters['leiden_resolution']))
    vis.pca_plot(True, 'leiden')

    # 2) UMAP
    vis.umap_alg()
    vis.umap_visualization(color='leiden', title='leiden')

    # 3) cell cycle inspection
    vis.cell_cycle_inspection(True)

    # 4) save dataset for CNV inference in R
    utils.save_adata_h5ad(adata, output_directory, filename, 'preprocessed')

    return adata, vis


def annotation_pipeline(adata, output_directory, filename):
    """
    Initial automated cell type annotation and malignant/non-malignant consensus call.
    """

    os.makedirs(output_directory + '/2.downstream_analysis', exist_ok=True)

    vis = visualization.Visualization(adata, output_directory + '/2.downstream_analysis', filename)

    print('-> Starting cell annotation: ')
    annotation = cell_annotation.Cell_Annotation(adata, output_directory + '/2.downstream_analysis', filename)

    # 1) marker-based annotation from differential expression
    annotation.differential_expr_genes()
    cluster_to_de = annotation.azimuth_annotation()

    # tumor marker expression
    annotation.tumor_marker_expression_sum()
    vis.umap_visualization(color='tumor_marker_expression', title='tumor_marker_expression')

    # 2) CellTypist
    cluster_to_celltypist = annotation.celltypist_prediction()

    # 3) PanglaoDB
    cluster_to_pangl_db = annotation.panglaodb_annotation()

    # 4) gene set enrichment
    cluster_to_gsea = annotation.gsea_annotation()

    # 5) visualize the four annotations before consensus
    vis.umap_visualization(color=['healthy_vs_tumor_DE', 'healthy_vs_tumor_celltypist',
                                  'healthy_vs_tumor_panglaodb', 'healthy_vs_tumor_gsea'],
                           title=['healthy_vs_tumor_DE', 'healthy_vs_tumor_celltypist',
                                  'healthy_vs_tumor_panglaodb', 'healthy_vs_tumor_gsea'])

    # 6) consensus annotation by majority vote
    annotation.consensus(cluster_to_de, cluster_to_celltypist, cluster_to_pangl_db, cluster_to_gsea)
    vis.umap_visualization(color='cluster_to_consensus_malignant', title='consensus malignant')
    vis.umap_visualization(color='cluster_to_consensus_all', title='consensus all')

    # 7) export cluster-level annotation for inferCNV
    adata.obs['malignant_vs_tumor_leiden'] = adata.obs['cluster_to_consensus_malignant'].astype(str) + '_' + \
                                             adata.obs['leiden'].astype(str)
    adata.obs['malignant_vs_tumor_leiden'].to_csv(output_directory + '/' + filename + '_malignant_annotation.csv',
                                                  sep='\t')

    # 8) save annotated dataset
    adata.obs = adata.obs.astype(str)
    utils.save_adata_h5ad(adata, output_directory, filename, 'malignant_vs_non_malignant')

    return adata


def main():
    parser = argparse.ArgumentParser(prog='Pipeline')
    parser.add_argument('-d', nargs=1, help='Path to MEX directory (filtered_feature_*) OR h5ad file', required=True)
    parser.add_argument('-j', nargs=1, help='.json file containing parameters')
    parser.add_argument('-n', nargs=1, help='Sample name used as prefix for output files. '
                                            'Default: name of the input file or directory')
    parser.add_argument('-o', nargs=1, help='Output directory. '
                                            'Default: <sample name>_data_analysis_<date> next to the input')
    parser.add_argument('-pr', help='Only Preprocessing', action='store_true')
    parser.add_argument('-dw', help='Only Downstream Analysis. File given must be h5ad processed', action='store_true')

    args = parser.parse_args()
    file_dir = os.path.abspath(args.d[0])
    only_dw = args.dw
    only_preprocessing = args.pr

    if args.j is not None:
        with open(args.j[0], 'r') as f:
            global parameters
            parameters = json.load(f)

    if not os.path.exists(file_dir):
        print(file_dir)
        parser.error("File or Directory containing data not found!")
    today = date.today().strftime("%b_%d_%Y")

    # sample name: taken from -n if given, otherwise from the input file or directory name
    if args.n is not None:
        filename = args.n[0]
    else:
        filename = os.path.basename(file_dir.rstrip('/'))
        filename = os.path.splitext(filename)[0]
        for suffix in ('_aligned', '_reanalyzed', '_preprocessed'):
            if filename.endswith(suffix):
                filename = filename[:-len(suffix)]

    # output directory: taken from -o if given, otherwise created next to the input
    if args.o is not None:
        output_directory = os.path.abspath(args.o[0])
    else:
        output_directory = os.path.join(os.path.dirname(file_dir),
                                        filename + '_data_analysis_' + today)
    os.makedirs(output_directory, exist_ok=True)

    # read data: var_names contains gene names, obs_names contains cell names
    print('Name of file: {}'.format(filename))
    print('-> Reading expression matrix ...')
    if os.path.isfile(file_dir):
        print('Reading Data File...')
        adata = sc.read_h5ad(file_dir)
    else:
        print('Reading Directory...')
        adata = sc.read_10x_mtx(file_dir, cache=True)
    nr_cells, nr_genes = adata.shape[0], adata.shape[1]
    print('   Number of cells (rows): {} || Number of genes (columns): {}'.format(adata.shape[0], adata.shape[1]))

    # STEP 1: preprocessing
    if not only_dw:
        adata = preprocessing(adata, output_directory, filename, nr_genes, nr_cells)

    # STEP 2: downstream analysis and annotation
    if not only_preprocessing:
        adata, vis = downstream_analysis(adata, output_directory, filename)
        adata.obs.to_csv(output_directory + '/metadata.csv')
        adata = annotation_pipeline(adata, output_directory, filename)

    print('Done. Results saved in directory ' + output_directory + '\n\n')


if __name__ == "__main__":
    print('Starting Analysis...')
    start = time.time()
    main()
    end = time.time()
    print('Time: {} min'.format((end - start) / 60))
