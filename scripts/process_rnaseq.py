"""
Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/process_rnaseq.py

Cleans up the downloaded RNAseq data and zero-one normalizes

Usage:
Run only once by run_pipeline.sh

Output:
Normalized PCL file with genes as rows and sample IDs as columns and MAD genes
if option is selected
"""

import os
import pandas as pd
from sklearn import preprocessing
from statsmodels.robust import scale
import argparse


def normalize_data(data, out_file, mad=False,
                   mad_file=os.path.join('tables', 'full_mad_genes.tsv'),
                   output=True,
                   method='minmax'):
    """
    Filters unidentified genes and normalizes each input gene expression matrix

    Arguments:
    :param data: pandas DataFrame genes as rows and sample IDs as columns
    :param out_file: the file name to write normalized matrix
    :param mad: boolean indicating if MAD genes should be output to file
    :param mad_fh: the file name to write mad genes
    :param method: the type of scaling to perform (defaults to minmax)

    Output:
    Writes normalized matrix (if output=True) and mad genes to file
    (if mad=True); returns the normalized matrix if output=False
    """

    # Drop all row names with unidentified gene
    data = data[-data.index.str.contains('?', regex=False)]

    # Sort data by gene name
    data = data.sort_index()

    # Zero-one normalize
    if method == 'minmax':
        min_max_scaler = preprocessing.MinMaxScaler()
        data_normalize = min_max_scaler.fit_transform(data.T)

    elif method == 'zscore':
        data_normalize = preprocessing.scale(data.T, axis=0)

    data_normalize = pd.DataFrame(data_normalize, index=data.columns,
                                  columns=data.index).T
    # Write to file
    if output:
        data_normalize.to_csv(out_file, sep='\t', header=True, index=True)
    else:
        return data_normalize

    # Write out MAD genes
    if mad:
        all_mad_genes = scale.mad(data_normalize, c=1, axis=1)
        all_mad_genes = pd.Series(all_mad_genes,
                                  index=data_normalize.index.values)
        all_mad_genes = all_mad_genes.sort_values(ascending=False)
        all_mad_genes.to_csv(mad_file, sep='\t', header=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tissue", dest="tissues",
                        help="tissues to normalize")
    parser.add_argument("-m", "--madgenes", dest="mad",
                        help="boolean to output mad genes", default=True)
    parser.add_argument("-c", "--concatenate", dest="concatenate",
                        help="boolean to normalize together", default=False)
    parser.add_argument("-s", "--scale", dest="scale_x",
                        help="method to scale x matrix", default='minmax')
    args = parser.parse_args()

    tissues = args.tissues.split(',')
    concatenate = args.concatenate
    mad = args.mad
    scale_x = args.scale_x

    if concatenate == 'all':
        pancan_file = args.tissues.replace(',', '-')
    elif concatenate:
        pancan_file = concatenate.replace(',', '-')

    # Load Data
    tissue_list = []
    for tissue in tissues:
        raw_file = os.path.join('data', 'X', 'raw',
                                '{}_RNAseq_X_matrix.tsv'.format(tissue))
        out_file = os.path.join('data', 'X', 'normalized',
                                '{}_X_normalized.tsv'.format(tissue))
        raw_df = pd.read_csv(raw_file, delimiter='\t', index_col=0)
        if concatenate:
            if concatenate == 'all':
                tissue_list.append(raw_df)
            elif tissue in concatenate.split(','):
                tissue_list.append(raw_df)
        else:
            normalize_data(raw_df, out_file, mad, method=scale_x)

    if concatenate:
        pancan_file = os.path.join('data', 'X', 'normalized',
                                   '{}_X_normalized.tsv'.format(pancan_file))
        raw_df = pd.concat(tissue_list, axis=1)
        normalize_data(raw_df, pancan_file)
