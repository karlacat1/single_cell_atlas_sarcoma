"""
Helper functions for saving AnnData objects.
"""

import os
import shutil
import subprocess
from scipy import io


def save_adata_h5ad(adata, output_directory, filename, file_ending):
    """
    Write an AnnData object to <output_directory>/<filename>_<file_ending>.h5ad
    """
    print('---> Saving ', file_ending, ' dataset...')
    adata.write(output_directory + '/' + filename + '_' + file_ending + '.h5ad')


def save_adata_for_r_script(adata_norm, output_dir):
    """
    Export an AnnData object as a gzipped 10x-style MEX triplet
    (barcodes.tsv, features.tsv, matrix.mtx) for downstream analysis in R.
    """
    print('Saving Matrix files ...')
    if os.path.exists(output_dir + '/matrix_files'):
        shutil.rmtree(output_dir + '/matrix_files')
    os.makedirs(output_dir + '/matrix_files')
    with open(output_dir + "/matrix_files/barcodes.tsv", 'w') as f:
        for item in adata_norm.obs_names:
            f.write(item + '\n')

    with open(output_dir + '/matrix_files/features.tsv', 'w') as f:
        for item in ['\t'.join([x, x, 'Gene Expression']) for x in adata_norm.var_names]:
            f.write(item + '\n')
    io.mmwrite(output_dir + '/matrix_files/matrix', adata_norm.X.T)
    gzip_cmd = 'gzip ' + output_dir + '/matrix_files/*'
    subprocess.run(gzip_cmd, shell=True)
