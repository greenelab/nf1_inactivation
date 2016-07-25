"""
(C) Gregory Way 2016
PanCancer NF1 Classifier
knn.py

Description:
Function to determine nearest neighbor of unlabeled data

Usage:
Import only by build_y.py
"""

import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier


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


def knn(X, mutation, n=2, gene=['NF1'], mut_filt=['Silent']):
    """
    Output classification for samples in X

    Arguments:
    :param X: pandas dataframe with genes as rows and samples as columns
    :param mut: pandas dataframe with mutation info
    :param n: the number of neighbors to subset
    :param gene: string indicating gene in mutation dataframe
    :param mut_filt: list of classes of mutations to filter

    Output:
    list of classifications for each sample
    0 - Not mutated
    1 - Mutated
    2 - Not assessed
    list of samples belonging to each assignment
    """

    # Process data
    y_class, zero_samp, one_samp, two_samp = get_y_class(X, mutation, gene,
                                                         mut_filt)
    x_clf = pd.concat([X[zero_samp], X[one_samp]], axis=1)
    x_test = X[two_samp]

    # generate y vector of form array([0, 0, ..., 0, 1, 1, ..., 1])
    y_ = np.array([0] * len(zero_samp) + [1] * len(one_samp))

    # Train a KNN classifier with known status
    neigh = KNeighborsClassifier(metric='euclidean', weights='distance',
                                 n_neighbors=int(n))
    neigh.fit(x_clf.T, y_)
    y_knn = neigh.predict(x_test.T)

    # Build the output
    samples = zero_samp + one_samp + two_samp
    labels = [0] * len(zero_samp) + [1] * len(one_samp) + list(y_knn)

    return pd.DataFrame([labels, samples]).T


def augment_knn(y, aug_x, n=2, max_iter=5):
    """
    Iterate through the data to augment the y values

    Arguments:
    :param Y: pandas DataFrame with [classification, sample ID] as columns
    :param X: pandas DataFrame with expression values
    :param n: number of neighbors
    :param max_iter: how many times to loop through the data
    """

    from random import shuffle

    y.columns = ['status', 'id']
    aug_y = y.copy()
    iteration = 0

    while iteration <= max_iter:
        sample_shuffle = list(range(aug_x.shape[1]))
        shuffle(sample_shuffle)

        y_old = aug_y.copy()
        for idx in sample_shuffle:
            samp = aug_x.columns[idx]

            # Get removed sample
            x_i = aug_x[samp]
            y_i = aug_y.loc[aug_y[aug_y['id'] == samp], :]
            y_idx = y_i.index.tolist()[0]
            y_i = y_i['status'].tolist()[0]

            # Because we know that mutations are gold standards
            # Perform augmented KNN for only class 0
            if y_i == 1:
                continue

            # Retrieve the full X and Y matrix without samp
            x_minus_i = aug_x.drop([samp], axis=1)
            y_minus_i = aug_y[aug_y['id'] != samp]
            y_minus_i = y_minus_i['status'].tolist()

            # Build classifications
            neigh = KNeighborsClassifier(metric='euclidean',
                                         weights='distance',
                                         n_neighbors=int(n))
            neigh.fit(x_minus_i.T, y_minus_i)
            y_hat = neigh.predict(x_i.T.reshape(1, -1))[0]

            # Assess the prediction
            # this will only change if y_i=0 and y_hat=1
            if y_i == y_hat:
                continue
            else:
                # Reassign the augmented y to y_hat
                aug_y.set_value(y_idx, 'status', y_hat)
        iteration += 1
        if y_old['status'].tolist() == aug_y['status'].tolist() and \
           iteration >= 2:
            break
    return aug_y
