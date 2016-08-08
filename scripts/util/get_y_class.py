"""
(C) Gregory Way 2016
PanCancer NF1 Classifier
scripts/util/get_y_class.py

Description:
Function to determine nearest neighbor of unlabeled data

Usage:
Import only by build_y.py
"""

import os
import pandas as pd


def get_y_class(mut, tis, gene=['NF1'], mut_filt=['Silent'],
                x_base='data/X/normalized/'):
    """
    Output classification for samples in X

    Arguments:
    :param mut: pandas dataframe with mutation info
    :param tis: the tissue to load X
    :param gene: string indicating gene in mutation dataframe
    :param mut_filt: list of classes of mutations to filter
    :param x_base: string of where the X matrices are stored

    Output:
    list of classifications for each sample
    0 - Not mutated
    1 - Mutated
    2 - Not assessed
    list of samples belonging to each assignment
    """
    # Read in X
    x_files = os.listdir(x_base)
    x_file = [x for x in x_files if x.split('_')[0] == tis]
    X = pd.read_csv(x_base + x_file[0], delimiter='\t', index_col=0)

    y_class = []
    zero_samp = []
    one_samp = []
    two_samp = []
    for samp in X.columns:
        # Subset mutation file to only consider that sample
        m = mut[mut['#sample'] == samp]

        # If the sample is not in the mutation file
        if m.shape[0] == 0:
            y_class.append(0)
            two_samp.append(samp)

        # If the sample is in the mutation file
        else:
            # Subset the subsetted mutation file to only consider the gene
            m = m.loc[m['gene'].isin(gene), :]

            # If the gene is not in the file, append 0
            if m.shape[0] == 0:
                y_class.append(0)
                zero_samp.append(samp)

            # Gene is in the file, but test if not filtered and append status
            else:

                if mut_filt is None:
                    filt = m

                else:
                    # Subset the mutation file again to silent mutations
                    filt = m.loc[m['effect'].isin(mut_filt), :]

                # If the shape is the same as the filtered, then append 0
                if m.shape[0] == filt.shape[0]:
                    y_class.append(0)
                    zero_samp.append(samp)

                # If not, then there is least one non-filtered mutation
                else:
                    y_class.append(1)
                    one_samp.append(samp)
    return y_class, zero_samp, one_samp, two_samp, X
