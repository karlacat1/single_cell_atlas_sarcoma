"""
Initial automated cell type annotation for a single sample.

Annotates Leiden clusters by four independent strategies and combines them into a
consensus call:

  1. marker_gene_overlap against a curated marker panel (`azimuth_annotation`)
  2. reference-based classification with CellTypist (`celltypist_prediction`)
  3. marker_gene_overlap against PanglaoDB signatures (`panglaodb_annotation`)
  4. gene set enrichment of ranked cluster markers against MSigDB (`gsea_annotation`)

Each strategy returns a per-cluster label and a coarse malignant / non-malignant
call; `consensus` derives a per-cluster consensus by majority vote.

This produces a starting point for manual curation. Final cell type and cell state
annotations used in the manuscript were assigned after cohort-wide integration.

External resource files are read from the directory given by the SARCOMA_RESOURCES
environment variable (default: ../resources relative to this file). See
resources/README.md for how to obtain them.
"""

import os

import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
import gseapy
import matplotlib.pyplot as plt
import seaborn as sns
import shutil

RESOURCES_DIR = os.environ.get(
    'SARCOMA_RESOURCES',
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'resources')
)

PANGLAO_FILE = os.path.join(RESOURCES_DIR, 'PanglaoDB_markers_27_Mar_2020.tsv')

# MSigDB v2023.1 human collections
GMT_CELLTYPE = os.path.join(RESOURCES_DIR, 'msigdb', 'c8.all.v2023.1.Hs.symbols.gmt')   # cell type signatures
GMT_CANCER_MODULES = os.path.join(RESOURCES_DIR, 'msigdb', 'c4.all.v2023.1.Hs.symbols.gmt')
GMT_ONCOGENIC = os.path.join(RESOURCES_DIR, 'msigdb', 'c6.all.v2023.1.Hs.symbols.gmt')
GMT_GENE_ONTOLOGY = os.path.join(RESOURCES_DIR, 'msigdb', 'c5.all.v2023.1.Hs.symbols.gmt')

GSEA_THREADS = 25

# Curated marker panel. Non-malignant entries are canonical lineage markers; the
# 'malignant' entry collects entity-associated markers across the sarcoma entities
# in the cohort.
immune_marker = {
    'Macrophages': ['CD163', 'CD14', 'CSF1R', 'ITGAL', 'ITGAM', 'CD68', 'FUT4', 'FCGR3A', 'CD33', 'CD11', 'CTSS',
                    'FCN1', 'NEAT1', 'LYZ', 'PSAP', 'S100A9', 'AIF1', 'MNDA', 'SERPINA1', 'TYROBP'],
    'CD4+ T cell': ['IL7R', 'MAL', 'LTB', 'CD4', 'LDHB', 'TPT1', 'TRAC', 'TMSB10', 'CD3D', 'CD3G'],
    'Monocyte': ['CTSS', 'FCN1', 'NEAT1', 'LYZ', 'PSAP', 'S100A9', 'AIF1', 'MNDA', 'SERPINA1', 'TYROBP'],
    'NK': ['NKG7', 'KLRD1', 'TYROBP', 'GNLY', 'FCER1G', 'PRF1', 'CD247', 'KLRF1', 'CST7', 'GZMB'],
    'CD8+ T cell': ['CD8B', 'CD8A', 'CD3D', 'TMSB10', 'HCST', 'CD3G', 'LINC02446', 'CTSW', 'CD3E', 'TRAC'],
    'B cell': ['CD79A', 'RALGPS2', 'CD79B', 'MS4A1', 'BANK1', 'CD74', 'TNFRSF13C', 'HLA-DQA1', 'IGHM', 'MEF2C'],
    'other T cells': ['CD3D', 'TRDC', 'GZMK', 'KLRB1', 'NKG7', 'TRGC2', 'CST7', 'LYAR', 'KLRG1', 'GZMA'],
    'DC': ['CD74', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'CCDC88A', 'HLA-DRA', 'HLA-DMA', 'CST3', 'HLA-DQB1',
           'HLA-DRB1'],
    'Fibroblasts': ['COL1A1', 'COL1A2', 'COL5A1', 'LOXL1', 'LUM', 'FBLN1', 'FBLN2', 'CD34', 'PDGFRA', 'IL1R1'],
    'Endothelial Cells': ['VWF', 'CD34', 'PECAM1', 'CD36', 'ENTPD1', 'CD44', 'ICAM1', 'CD47', 'ITGB3', 'CD80'],
    'Epithelial cells': ['MUC1', 'ALCAM', 'ANPEP', 'ACE', 'RNPEP', 'MUC16', 'CD46', 'CEACAM1', 'KRT7', 'ITGB4'],
    'malignant': ['ALK', 'APC', 'CD99', 'NCAM1', 'DES', 'VIM', 'BCL2', 'B2M', 'KIT', 'SLC17A8', 'MAPT-AS1', 'POSTN',
                  'ADAM12', 'HMCN1', 'GRIA4', 'D21S2088E', 'CCL1', 'HTN3', 'GPR182', 'MADCAM1', 'FLT4', 'SCGB3A1',
                  'NTS', 'MMP12', 'CTLA4', 'CST5', 'ARHGAP36', 'CDH10', 'TMEM215', 'IGSF1', 'ZCCHC12', 'CST2',
                  'GAL3ST3', 'UNC5D', 'FOXL2NB', 'PCGEM1', 'SLC17A8', 'SCNN1G', 'C6orf15', 'PGLYRP2', 'CALCA',
                  'C1QL2', 'LIPI', 'MAPT-AS1', 'NKX2-2', 'SERPINB11', 'LAMP5', 'ADAM12', 'POSTN', 'ADAMTS12',
                  'XPNPEP2', 'PLA2G2A', 'F2RL1', 'ZNF469', 'HMCN1', 'GJB6', 'GFAP', 'GAP43', 'ANGPTL7', 'SOX10',
                  'GRIA4', 'SERPINA5', 'ITGB8', 'S100B', 'FABP7', 'C5orf58', 'SPATA22', 'ZIC1', 'FGF19', 'GPR126',
                  'ZIC4', 'FAM163B', 'GFRA2', 'DLX6-AS1', 'HOXC11', 'TOX', 'TP53', 'PAX3', 'FOXO1', 'FGFR4']
}

# PanglaoDB cell types considered for annotation
PANGLAO_CELLTYPES = ['B cells', 'B cells memory', 'B cells naive', 'NK cells', 'Endothelial cells', 'Macrophages',
                     'Epithelial cells', 'Dendritic cells', 'Fibroblasts', 'Monocytes', 'T cells',
                     'Gamma delta T cells', 'Natural killer T cells', 'T cells naive', 'T cytotoxic cells',
                     'T follicular helper cells', 'T helper cells', 'T memory cells', 'T regulatory cells']


class Cell_Annotation:
    def __init__(self, adata, output_directory, filename):
        """
        :param adata: preprocessed and clustered AnnData object for one sample
        :param output_directory: directory for annotation plots and tables
        :param filename: sample name, used as a prefix for output files
        """
        self.adata_norm = adata
        self.output_directory = output_directory
        self.filename = filename
        self.clusters = np.unique(adata.obs['leiden'].astype('int'))

    def differential_expr_genes(self):
        """
        Rank marker genes per Leiden cluster with a Wilcoxon rank-sum test and
        write the top 15 markers per cluster.
        """
        print('---> Computing DE genes ')

        sc.tl.rank_genes_groups(self.adata_norm, groupby='leiden', method='wilcoxon')
        sc.pl.rank_genes_groups(self.adata_norm, n_genes=15, show=False, sharey=False)
        plt.savefig(self.output_directory + '/' + self.filename + '_DE_genes_wilcoxon.png', bbox_inches='tight')
        plt.close()
        self.adata_norm.wilcoxon = self.adata_norm.uns['rank_genes_groups']

        stat_test_top = pd.DataFrame(self.adata_norm.wilcoxon['names'])[:15]
        stat_test_top.to_csv(self.output_directory + '/' + self.filename + '_DE_top_marker.csv', sep='\t')

    def azimuth_annotation(self):
        """
        Assign each cluster to the cell type in the curated marker panel with the
        highest marker overlap. Clusters with no overlap are called malignant;
        clusters tying between a malignant and a non-malignant label are unclear.
        """
        print('---> Annotating based on markers')

        marker_matches = sc.tl.marker_gene_overlap(self.adata_norm, immune_marker, adj_pval_threshold=0.05)
        marker_matches = marker_matches[marker_matches.sum(axis=1) > 0]

        fig, ax = plt.subplots(figsize=(16, 8))
        sns.heatmap(marker_matches, linewidth=0.2)
        plt.savefig(self.output_directory + '/' + self.filename + '_DE_celltype_heatmap.png',
                    dpi=100, bbox_inches='tight')
        plt.close()

        cluster_to_de = dict(marker_matches.idxmax())

        for cluster in marker_matches.columns:
            max_overlap = np.max(marker_matches[cluster])
            # no overlap with any cell type -> possible malignant
            if max_overlap == 0:
                cluster_to_de[cluster] = 'malignant'
                continue

            max_celltype = marker_matches.index[(marker_matches[cluster] == max_overlap)]
            # tie between a malignant and a non-malignant label -> unclear
            if len(max_celltype) > 1:
                if 'malignant' in max_celltype:
                    cluster_to_de[cluster] = 'unclear'

        self.adata_norm.obs['cluster_to_de'] = self.adata_norm.obs['leiden'].map(cluster_to_de)

        cluster_contains_all = cluster_to_de.copy()
        # collapse all non-malignant labels for the coarse call
        for i in cluster_to_de.keys():
            if cluster_to_de[i] != 'malignant' and cluster_to_de[i] != 'unclear':
                cluster_to_de[i] = 'non malignant'

        self.adata_norm.obs['healthy_vs_tumor_DE'] = self.adata_norm.obs['leiden'].map(cluster_to_de)
        nr_non_malignant = np.sum(self.adata_norm.obs['healthy_vs_tumor_DE'] == 'non malignant')
        nr_malignant = np.sum(self.adata_norm.obs['healthy_vs_tumor_DE'] == 'malignant')
        unclear = np.sum(self.adata_norm.obs['healthy_vs_tumor_DE'] == 'unclear')

        print('% of non malignant cells:', nr_non_malignant / self.adata_norm.shape[0] * 100)
        print('% of malignant cells:', nr_malignant / self.adata_norm.shape[0] * 100)
        print('% of unclear cells:', unclear / self.adata_norm.shape[0] * 100)

        return cluster_contains_all

    def tumor_marker_expression_sum(self):
        """
        Summed expression of the entity-associated marker panel per cell, computed
        on the full gene set stored in .raw.
        """
        tumor_marker = np.unique(immune_marker['malignant'])
        adata_raw = self.adata_norm.raw.to_adata()
        markers_in_data = [i in adata_raw.var.index.values for i in tumor_marker]
        markers_in_data = tumor_marker[markers_in_data]

        adata_tumor = adata_raw[:, markers_in_data]
        t_sum = np.sum(pd.DataFrame(adata_tumor.X.toarray()), axis=1)
        self.adata_norm.obs['tumor_marker_expression'] = t_sum.values

    def celltypist_prediction(self, model='Immune_All_High.pkl'):
        """
        Classify cells with CellTypist using majority voting over Leiden clusters.
        Cells predicted with a confidence score below 0.99 are treated as
        malignant, since the reference model covers immune cell types only.
        """
        print('---> Annotating using Celltypist')

        predictions = celltypist.annotate(self.adata_norm, model=model, majority_voting=True)
        self.adata_norm = predictions.to_adata()

        self.adata_norm.obs['cells'] = self.adata_norm.obs.index.values

        self.adata_norm.obs['highest_confidence'] = self.adata_norm.obs['majority_voting'].astype('str')
        self.adata_norm.obs.loc[
            self.adata_norm.obs['conf_score'].astype('double') < 0.99, "highest_confidence"] = 'malignant'

        celltypes = np.unique(self.adata_norm.obs["highest_confidence"])
        cluster_to_celltype_conf = {}
        celltype_distribution_conf = pd.DataFrame(0.0, columns=self.clusters, index=celltypes)
        for name, group in self.adata_norm.obs.groupby(by='leiden'):
            unique, counts = np.unique(group['highest_confidence'], return_counts=True)
            cluster_to_celltype_conf[name] = unique[np.argmax(counts)]
            celltype_distribution_conf.loc[unique, int(name)] = counts / np.sum(counts) * 100

        celltype_distribution_conf.T.plot(kind='bar', stacked=True, figsize=(14, 8))
        plt.title('Cell type distribution pro cluster:')
        plt.legend(bbox_to_anchor=(1.0, 1.0))
        plt.savefig(self.output_directory + '/' + self.filename + '_confidence_score.png', bbox_inches='tight')
        plt.close()

        self.adata_norm.obs['cluster_to_celltypist'] = self.adata_norm.obs['leiden'].astype(str).replace(cluster_to_celltype_conf)
        self.adata_norm.obs['healthy_vs_tumor_celltypist'] = self.adata_norm.obs['cluster_to_celltypist'].astype(
            'string')
        self.adata_norm.obs.loc[
            self.adata_norm.obs['healthy_vs_tumor_celltypist'] != "malignant",
            "healthy_vs_tumor_celltypist"] = 'non malignant'

        nr_non_malignant = np.sum(self.adata_norm.obs['healthy_vs_tumor_celltypist'] == 'non malignant')
        nr_malignant = np.sum(self.adata_norm.obs['healthy_vs_tumor_celltypist'] != 'non malignant')

        print('Nr of malignant clusters:', np.sum(np.asarray(list(cluster_to_celltype_conf.values())) == 'malignant'))
        print('Nr of non malignant clusters:',
              np.sum(np.asarray(list(cluster_to_celltype_conf.values())) != 'malignant'))
        print('Nr of non malignant cells:', nr_non_malignant)
        print('Nr of malignant cells:', nr_malignant)

        return cluster_to_celltype_conf

    def panglaodb_annotation(self):
        """
        Assign each cluster to the PanglaoDB cell type signature with the highest
        marker overlap, using the same malignant / unclear rules as
        `azimuth_annotation`.
        """
        print('---> Annotating using PanglaoDB')
        pangla_db = pd.read_csv(PANGLAO_FILE, sep='\t')

        gene_marker_dict = {}
        for current_celltype in PANGLAO_CELLTYPES:
            gene_marker = pangla_db[pangla_db['cell type'] == current_celltype]['official gene symbol']
            gene_marker_dict[current_celltype] = gene_marker
        gene_marker_dict['malignant'] = pd.Series(immune_marker['malignant'])

        marker_matches_pangl = sc.tl.marker_gene_overlap(self.adata_norm, gene_marker_dict)
        marker_matches_pangl = marker_matches_pangl[marker_matches_pangl.sum(axis=1) != 0]

        cluster_to_pangl_db = dict(marker_matches_pangl.idxmax())

        for cluster in marker_matches_pangl.columns:
            max_overlap = np.max(marker_matches_pangl[cluster])

            if max_overlap == 0:
                print(cluster)
                cluster_to_pangl_db[cluster] = 'malignant'
                continue
            max_celltype = marker_matches_pangl.index[(marker_matches_pangl[cluster] == max_overlap)]
            if len(max_celltype) > 1:
                if 'malignant' in max_celltype:
                    cluster_to_pangl_db[cluster] = 'unclear'

        self.adata_norm.obs['cluster_to_panglaodb'] = self.adata_norm.obs['leiden'].astype(str).replace(cluster_to_pangl_db)
        cluster_to_complete = dict(cluster_to_pangl_db)

        for i in cluster_to_pangl_db.keys():
            if cluster_to_pangl_db[i] != 'malignant' and cluster_to_pangl_db[i] != 'unclear':
                cluster_to_pangl_db[i] = 'non malignant'

        self.adata_norm.obs['healthy_vs_tumor_panglaodb'] = self.adata_norm.obs['leiden'].astype(str).replace(cluster_to_pangl_db)
        nr_non_malignant = np.sum(self.adata_norm.obs['healthy_vs_tumor_panglaodb'] == 'non malignant')
        nr_malignant = np.sum(self.adata_norm.obs['healthy_vs_tumor_panglaodb'] != 'non malignant')

        print('Nr of malignant clusters:', np.sum(np.asarray(list(cluster_to_pangl_db.values())) == 'malignant'))
        print('Nr of non malignant clusters:', np.sum(np.asarray(list(cluster_to_pangl_db.values())) != 'malignant'))
        print('Nr of non malignant cells:', nr_non_malignant)
        print('Nr of malignant cells:', nr_malignant)

        return cluster_to_complete

    def gsea_annotation(self):
        """
        Rank each cluster's markers by log fold change and run preranked GSEA
        against four MSigDB collections. Clusters whose top gene sets come mainly
        from the cancer module (c4) or oncogenic signature (c6) collections are
        called malignant; otherwise non-malignant.
        """
        if os.path.exists(self.output_directory + '/gsea'):
            shutil.rmtree(self.output_directory + '/gsea')
        os.makedirs(self.output_directory + '/gsea')

        results = {}
        for cl in np.unique(self.adata_norm.obs['leiden'].astype('str')):
            results[cl] = 'non malignant'

        for cluster in np.unique(self.adata_norm.obs['leiden']):
            print('Processing of cluster ', cluster)

            # rank genes by log fold change, keeping genes detected in >30 cells
            gene_rank = sc.get.rank_genes_groups_df(self.adata_norm, group=cluster)[['names', 'logfoldchanges']]
            gene_rank.sort_values(by=['logfoldchanges'], inplace=True, ascending=False)
            self.adata_norm.var.n_cells_by_counts = self.adata_norm.var.n_cells_by_counts.astype('float64')
            gene_rank = gene_rank[
                gene_rank['names'].isin(self.adata_norm.var_names[self.adata_norm.var.n_cells_by_counts > 30])]
            print('Gene ranks done.')

            print('Applying GSEA.')
            res = gseapy.prerank(rnk=gene_rank,
                                 gene_sets=[GMT_CELLTYPE, GMT_CANCER_MODULES, GMT_ONCOGENIC, GMT_GENE_ONTOLOGY],
                                 threads=GSEA_THREADS, permutation_num=1000, seed=6)
            terms = res.res2d.Term
            plot_range = 5
            if len(terms) == 0:
                print('No overlaps for cluster ', cluster)
                results[cluster] = 'malignant'
                continue
            nr_of_top_genes = 5
            if len(terms) < nr_of_top_genes:
                print('Low number of overlaps.')
                nr_of_top_genes = len(terms)
                if len(terms) < 3:
                    plot_range = len(terms)

            print('Computing gene set with high overlap.')
            # which collection do the top gene sets come from
            first_ten = terms.apply(lambda x: x.split('.')[0])[:nr_of_top_genes]
            unique, counts = np.unique(first_ten, return_counts=True)
            if pd.DataFrame(counts, index=unique).idxmax()[0] == 'c4' \
                    or pd.DataFrame(counts, index=unique).idxmax()[0] == 'c6':
                print(cluster, ' -> malignant')
                results[cluster] = 'malignant'
            else:
                print(cluster, ' -> possibly non malignant')
                results[cluster] = 'non malignant'

            print('Plotting...')
            for i in range(plot_range):
                gseapy.gseaplot(rank_metric=res.ranking, term=terms[i], **res.results[terms[i]],
                                ofname=self.output_directory + '/gsea/' + 'gsea_' + cluster + '_' + str(i) + '.png')

        self.adata_norm.obs['healthy_vs_tumor_gsea'] = self.adata_norm.obs['leiden'].map(results).astype('category')
        print('Nr of malignant clusters:', np.sum(np.asarray(list(results.values())) == 'malignant'))
        print('Nr of non malignant clusters:', np.sum(np.asarray(list(results.values())) != 'malignant'))

        return results

    def consensus(self, cluster_to_de, cluster_to_celltypist, cluster_to_panglaodb, cluster_to_gsea):
        """
        Combine the four per-cluster annotations by majority vote, once at the
        level of cell type labels and once at the level of the coarse
        malignant / non-malignant call.
        """

        def consensus_celltype(x):
            unique, counts = np.unique(x, return_counts=True)
            x_mod = np.asarray(['non malignant' if (i != 'unclear' and i != 'malignant') else i for i in x])
            unique_mod, counts_mod = np.unique(x_mod, return_counts=True)

            if np.max(counts) >= 3:
                return unique[np.argmax(counts)]
            elif np.max(counts_mod) >= 3:
                if np.max(counts) == 2:
                    return unique[np.argmax(counts)]
                else:
                    return unique_mod[np.argmax(counts_mod)]
            elif np.max(counts) == 2 and len(unique_mod) == 3:
                return unique[np.argmax(counts)]
            else:
                return 'unclear'

        consensus = pd.DataFrame(
            [cluster_to_de.values(), cluster_to_celltypist.values(), cluster_to_panglaodb.values(),
             cluster_to_gsea.values()],
            index=['DE', 'Celltypist', 'PanglaoDB', 'GSEA']).T

        # harmonize label names across the four sources
        consensus = consensus.replace('possible malignant', 'malignant')
        consensus = consensus.replace('T memory cells', 'T cells')
        consensus = consensus.replace('Gamma delta T cells', 'T cells')
        consensus = consensus.replace('CD4+ T cell', 'T cells')
        consensus = consensus.replace('Naive B cells', 'B cells')
        consensus = consensus.replace('Dendritic cells', 'DC')
        consensus = consensus.replace('Endothelial Cells', 'Endothelial cells')
        consensus = consensus.replace('Tcm/Naive helper T cells', 'T cells')
        consensus['results'] = consensus.apply(lambda x: consensus_celltype(x), axis=1)
        consensus.to_csv(self.output_directory + '/' + self.filename + '_consensus_result_complete.csv', sep='\t')
        cluster_to_consensus = dict(zip(consensus.index.astype('str'), consensus['results']))
        self.adata_norm.obs['cluster_to_consensus_all'] = self.adata_norm.obs['leiden'].astype(str).replace(cluster_to_consensus)

        consensus = pd.DataFrame(
            [cluster_to_de.values(), cluster_to_celltypist.values(), cluster_to_panglaodb.values(),
             cluster_to_gsea.values()],
            index=['DE', 'Celltypist', 'PanglaoDB', 'GSEA']).T
        consensus[np.logical_and(consensus != 'malignant', consensus != 'unclear')] = 'non malignant'
        consensus['results'] = consensus.apply(lambda x: consensus_celltype(x), axis=1)
        consensus.to_csv(self.output_directory + '/' + self.filename + '_consensus_result.csv', sep='\t')

        cluster_to_consensus = dict(zip(consensus.index.astype('str'), consensus['results']))
        self.adata_norm.obs['cluster_to_consensus_malignant'] = self.adata_norm.obs['leiden'].astype(str).replace(
            cluster_to_consensus)

        print('Nr of cells: ', self.adata_norm.shape[0])
        print('Nr of malignant cells: ', np.sum(self.adata_norm.obs['cluster_to_consensus_malignant'] == 'malignant'))
        print('Nr of non malignant cells: ',
              np.sum(self.adata_norm.obs['cluster_to_consensus_malignant'] == 'non malignant'))
        print('Nr of unclear cells: ', np.sum(self.adata_norm.obs['cluster_to_consensus_malignant'] == 'unclear'))
        pd.DataFrame({
            'Nr_cells': [self.adata_norm.shape[0]],
            'malignant': [np.sum(self.adata_norm.obs['cluster_to_consensus_malignant'] == 'malignant')],
            'non_malignant': [np.sum(self.adata_norm.obs['cluster_to_consensus_malignant'] == 'non malignant')],
            'unclear': [np.sum(self.adata_norm.obs['cluster_to_consensus_malignant'] == 'unclear')]
        }).to_csv(self.output_directory + '/' + self.filename + '_count_result.csv', sep='\t')
