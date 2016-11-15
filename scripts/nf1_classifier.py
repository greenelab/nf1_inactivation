"""
Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/nf1_classifier.py

Performs cross validation to select optimal hyperparameters.

Usage:
The script is highly modular and can be used in a variety of ways depending on
the flags. The script is run by run_pipeline.sh with a specific set of
parameters, but it can be built upon.

    With required flags:

    --gene             The gene to determine mutation status for (e.g. 'NF1')
    --tissue           Tissue acronym (e.g. GBM for glioblastoma)
    --ouput-file       The location to save the cross validation results
    --classifier       Classifier to use (can be either 'elasticnet' or 'SVM')
    --y-origin-fh      The location where the Y matrix is stored
    --hyperparameters  The hyperparameters to test during cross validation

    And the following optional flags:

    --alt-params       A second set of hyperparameters to test
    --folds            The number of cross validation folds to compute
    --holdout          Option to extract a holdout set (defaults to False)
    --fraction-hold    What percentage of data to holdout (defaults to 10%)
    --min-samples      Minimum samples with mutation in tissue (defaults to 5)
    --iterations       Number training/testing to perform (defaults to 50)
    --num_mad_genes    Number of MAD genes to subset the X (defaults to 8000)
    --x-origin-fh      The location where the X matrix is stored
    --plot-roc         Should an ROC curve be output
    --validation-sub   Subset based on validation genes (defaults to False)
"""

import os
import random
import argparse
import pandas as pd
import numpy as np
from subprocess import call

from sklearn.linear_model import SGDClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, roc_curve

from util.cancer_cv import cancer_cv

# Set seed
random.seed(123)

####################################
# Load Functions
####################################


def generate_y_dict(sample_dictionary, gene, tissue, y_base='data/Y/'):
    '''
    Generates a dictionary of tissue specific y vectors

    Arguments:
    :param sample_dictionary: tissues as keys and values are sample IDs
    :param gene: the focus of the y vector of labels
    :param tissue: a list of tissues to use to generate y_dict
    :param y_base: where the Y matrices are stored

    Output:
    Dictionary with tissues as keys and mutation status DataFrame as values
    '''

    ####################################
    # Load X and y
    ####################################
    y_files = os.listdir(y_base)
    y_files = [y for y in y_files if y != '.empty']

    y_dict = {}
    for tis in sample_dictionary.keys():
        if tis in tissue:
            y_file = [y for y in y_files if y.split('_')[1] == tis]
            y_file = [y for y in y_file if y.split('_')[2][:-4] == gene]
            y_file = y_base + y_file[0]

            # Read in the file and generate dataframe
            y = pd.read_csv(y_file, delimiter='\t')

            if y.shape[1] == 1:
                sample_names = sample_dictionary[tis]
                y['Samples'] = sample_names

            y['Tissue'] = np.array([tis] * y.shape[0])
            y.columns = ['Status', 'Samples', 'Tissue']
            y.index = y['Samples'].tolist()
            y_dict[tis] = y

    return y_dict


def train_test_error(clf, train_x, train_y, test_x, test_y):
    """
    Return training and testing errors according to input classifier

    Arguments:
    :param clf: classifier sklearn object
    :param train_x: gene expression matrix
    :param train_y: list of labels
    :param test_x: gene expression matrix
    :param test_y: list of labels

    Output:
    Returns training and testing auroc
    """
    model = clf.fit(train_x, train_y)
    pred_y = model.predict(train_x)
    train_err = roc_auc_score(train_y, pred_y, average='weighted')
    pred_y = model.predict(test_x)
    test_err = roc_auc_score(test_y, pred_y, average='weighted')
    return train_err, test_err

####################################
# Perform Analysis
####################################
if __name__ == '__main__':
    ####################################
    # Load Command Arguments
    ####################################
    parser = argparse.ArgumentParser()  # Load command line options

    # Required parameters
    parser.add_argument("-g", "--gene", dest="gene",
                        help="which gene to prepare cross validation")
    parser.add_argument("-t", "--tissue", dest="tissue",
                        help="which X matrix to load")
    parser.add_argument("-o", "--output-file", dest="out_fh",
                        help="location to write the results to")
    parser.add_argument("-c", "--classifier", dest="classifier",
                        help="the type of classifier to use")
    parser.add_argument("-e", "--y-origin-fh", dest="y_origin_fh",
                        help="Where the y matrices are pulled from",
                        default='data/Y/')
    parser.add_argument("-y", "--hyperparameters", dest="hyperparameters",
                        help="hyperparameters to optimize in cross validation")

    # Optional parameters
    parser.add_argument("-x", "--alt-params", dest="alt_params",
                        help="additional set of hyperparameters to optimize in"
                             "cross validation", default=0.1)
    parser.add_argument("-f", "--num_folds", dest="num_folds", default=5,
                        help="the number of cross validation fold intervals")
    parser.add_argument("-z", "--holdout", dest="holdout", default=False,
                        help="boolean of extracting holdout set or not")
    parser.add_argument("-p", "--fraction-holdout", dest="fraction_holdout",
                        help="Fraction of samples to be heldout", default=0.1)
    parser.add_argument("-m", "--min-samples", dest="min_samples",
                        help="the minimum number of samples with a mutation to"
                        "consider the tissue", default=5)
    parser.add_argument("-i", "--iterations", dest="iterations",
                        help="number of random seed shuffles", default=50)
    parser.add_argument("-n", "--num_mad_genes", dest="num_mad_genes",
                        help="number of mad genes to subset", default=8000)
    parser.add_argument("-d", "--x-origin-fh", dest="x_origin_fh",
                        help="Where the x matrices are pulled from",
                        default='data/X/normalized/')
    parser.add_argument("-r", "--plot-roc", dest="roc",
                        help="decide if an ROC plot should be output",
                        action='store_true')
    parser.add_argument("-v", "--validation-sub", dest='validation_sub',
                        action='store_true')
    parser.add_argument("-w", "--coef-weights", dest='coefficients',
                        action='store_true')
    args = parser.parse_args()

    ####################################
    # Load Constants
    ####################################
    # Required parameters
    GENE = args.gene
    TISSUE = args.tissue.replace(',', '-')
    OUT_FH = args.out_fh
    CLASSIFIER = args.classifier
    Y_BASE = args.y_origin_fh
    C = args.hyperparameters.split(',')

    # Optional parameters
    L = args.alt_params.split(',')
    NUM_FOLDS = args.num_folds
    HOLD = args.holdout
    FRACTION_HOLDOUT = args.fraction_holdout
    MIN_SAMPLES = args.min_samples
    K = args.iterations
    RANDOM_SEEDS = random.sample(range(1, 1000000), int(K))
    NUM_MAD_GENES = args.num_mad_genes
    X_BASE = args.x_origin_fh
    ROC = args.roc
    VALIDATION_SUBSET = args.validation_sub

    # Generate filenames to save output plots
    FIG_BASE = OUT_FH.split('/')[1].split('.')[0]
    if ROC:
        OUTPUT_FIGURE = 'figures/' + TISSUE + '_' + CLASSIFIER + '_' + \
                         FIG_BASE + 'alpha_' + str(C[0]) + '_l1ratio_' + \
                         str(L[0]) + '_roc.png'
    else:
        OUTPUT_FIGURE = 'figures/' + TISSUE + '_' + CLASSIFIER + '_' + \
                        FIG_BASE + '_cv.png'

    # Are we interested in the weights/coefficients/genes?
    GET_COEF = args.coefficients
    WEIGHT_FH = 'results/' + TISSUE + '_' + CLASSIFIER + '_' + FIG_BASE + \
                '_weights.tsv'

    ####################################
    # Load and Process Data
    ####################################
    # RNAseq dictionary with sample IDs in order of columns
    sample_dictionary = pd.read_pickle('output/RNAseq_sample_dictionary.p')

    # Generate Y Dictionary
    y_dict = generate_y_dict(sample_dictionary, GENE, TISSUE, Y_BASE)

    # Load X matrix
    if X_BASE == 'data/X/normalized/':
        x_fh = X_BASE + TISSUE + '_X_normalized.tsv'
    else:
        x_fh = X_BASE + TISSUE + '.tsv'

    X = pd.read_csv(x_fh, delimiter='\t', index_col=0)

    # subset using MAD genes
    mad = pd.Series.from_csv('tables/full_mad_genes.tsv', sep='\t')
    mad = mad[:NUM_MAD_GENES]
    X = X.ix[mad.index]
    X = X.dropna()

    # Build a DataFrame of zeros for storing weights if option is selected
    empty_zero = np.zeros(X.shape[0])
    weight_coef = pd.DataFrame(empty_zero, index=X.index, columns=['weight'])

    if VALIDATION_SUBSET:
        V = pd.read_csv('data/validation/normalized/validation_set.tsv',
                        delimiter='\t', index_col=0)
        V = V.ix[mad.index]
        V = V.dropna()
        X = X.ix[V.index]

    # roc output is iteratively made if ROC curve is to be output
    if ROC:
        col_names = ['tpr', 'fpr', 'auc', 'type']
        roc_output = pd.DataFrame(columns=col_names)

    ####################################
    # Generate Cross Validation Parameters
    ####################################
    # Do this for each random seed
    cv_accuracy = []
    for seed in RANDOM_SEEDS:
        print(seed)
        # Build a cross validation partition object
        pancan_cv = cancer_cv(y_dict, NUM_FOLDS, HOLD, MIN_SAMPLES, seed)

        # Determine fraction to holdout and assign folds for y and x matrices
        if HOLD:
            pancan_cv.holdout_samples(FRACTION_HOLDOUT)

        # Ensures equal partitioning into training and testing datasets
        pancan_cv.assign_folds()
        pancan_cv.assign_x_matrix_folds(X)

        # Perform NUM_FOLDS cross validation on the total dataset
        for fold in range(0, NUM_FOLDS):

            # Split according to the testing fold
            x_train, x_test, y_train, y_test = pancan_cv.split_train_test(fold)

            for param in C:
                if CLASSIFIER == 'SVM':
                    clf = SVC(C=float(param), class_weight='balanced')
                    tr_err, te_err = train_test_error(clf, x_train, y_train,
                                                      x_test, y_test)
                    cv_accuracy.append(['train', tr_err, param, TISSUE, seed])
                    cv_accuracy.append(['test', te_err, param, TISSUE, seed])

                elif CLASSIFIER == 'elasticnet':
                    for param_b in L:
                        clf = SGDClassifier(loss='log',
                                            penalty='elasticnet',
                                            alpha=float(param),
                                            l1_ratio=float(param_b),
                                            class_weight='balanced',
                                            random_state=seed)
                        if ROC:
                            clf.fit(x_train, y_train)
                            y_score = clf.decision_function(x_test)
                            y_score_train = clf.decision_function(x_train)

                            fpr, tpr, _ = roc_curve(y_test, y_score)
                            fpr_t, tpr_t, _ = roc_curve(y_train, y_score_train)

                            test_auc = roc_auc_score(y_test, y_score,
                                                     average='weighted')
                            train_auc = roc_auc_score(y_train, y_score_train,
                                                      average='weighted')

                            temp_test = pd.DataFrame([tpr, fpr,
                                                     [test_auc] * len(tpr),
                                                     ['test'] * len(tpr)],
                                                     index=col_names).T
                            temp_train = pd.DataFrame([tpr_t, fpr_t,
                                                      [train_auc] * len(tpr_t),
                                                      ['train'] * len(tpr_t)],
                                                      index=col_names).T

                            full = temp_test.append(temp_train,
                                                    ignore_index=True)
                            roc_output = roc_output.append(full,
                                                           ignore_index=True)
                            if GET_COEF:
                                p_add = pd.DataFrame(clf.coef_.T,
                                                     index=x_test.columns,
                                                     columns=['weight'])
                                weight_coef = weight_coef.add(p_add,
                                                              fill_value=0)

                        else:

                            train_err, test_err = train_test_error(clf,
                                                                   x_train,
                                                                   y_train,
                                                                   x_test,
                                                                   y_test)
                            cv_accuracy.append(['train', train_err, param,
                                                param_b, TISSUE, seed])
                            cv_accuracy.append(['test', test_err, param,
                                                param_b, TISSUE, seed])

    ####################################
    # Write results to file
    ####################################
    if ROC:
        roc_output.to_csv(OUT_FH, sep='\t', index=False)
    else:
        cv_accuracy = pd.DataFrame(cv_accuracy)
        if CLASSIFIER == 'SVM':
            cv_accuracy.columns = ['class', 'error', 'param', 'tissue', 'seed']
        elif CLASSIFIER == 'elasticnet':
            cv_accuracy.columns = ['class', 'error', 'alpha', 'l1_ratio',
                                   'tissue', 'seed']
        cv_accuracy.to_csv(OUT_FH, index=False, sep='\t')

    if GET_COEF:
        weight_coef['abs'] = weight_coef['weight'].abs()
        sorted_weight = weight_coef.sort_values('abs', ascending=False)
        sorted_weight = sorted_weight['weight']
        sorted_weight.to_csv(WEIGHT_FH, index=True, sep='\t')

        command = 'R --no-save --args ' + WEIGHT_FH + ' ' + \
                  'figures/weight_genes.png < scripts/viz/viz_geneweights.r'
        call(command, shell=True)

    ####################################
    # Plot
    ####################################
    if ROC:
        command = 'R --no-save --args ' + OUT_FH + ' ' + OUTPUT_FIGURE + \
                  ' < ' + 'scripts/viz/viz_roc.r'
    else:
        command = 'R --no-save --args ' + OUT_FH + ' ' + OUTPUT_FIGURE + \
                  ' ' + CLASSIFIER + ' < ' + 'scripts/viz/viz_cv.r'
    call(command, shell=True)
