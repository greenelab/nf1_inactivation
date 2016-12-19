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
    --ouput_file       The location to save the cross validation results
    --classifier       Classifier to use (can be either 'elasticnet' or 'SVM')
    --y_origin_fh      The location where the Y matrix is stored
    --hyperparameters  The hyperparameters to test during cross validation

    And the following optional flags:

    --alt_params       A second set of hyperparameters to test
    --folds            The number of cross validation folds to compute
    --holdout          Option to extract a holdout set (defaults to False)
    --fraction_hold    What percentage of data to holdout (defaults to 10%)
    --min_samples      Minimum samples with mutation in tissue (defaults to 5)
    --iterations       Number training/testing to perform (defaults to 50)
    --num_mad_genes    Number of MAD genes to subset the X (defaults to 8000)
    --x_origin_fh      The location where the X matrix is stored
    --plot_roc         Should an ROC curve be output
    --validation_sub   Subset based on validation genes (defaults to False)
    --effect_size      Calculate effect size of the cross validation intervals
"""

import os
import random
import argparse
import pandas as pd
import numpy as np
from subprocess import call
import matplotlib.pyplot as plt
import seaborn as sns

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


def get_cohens_d(score_df):
    """
    Find the effect size of decision functions for the input scoring dataframe
    Modified from:
    https://github.com/AllenDowney/CompStats/blob/master/effect_size.ipynb

    Arguments:
    :param score_df: pandas DataFrame with decision function scores and status

    Output:
    Returns the Cohen's d effect size estimate statistic for two status groups
    """
    pos_samples = np.array(score_df[score_df['status'] == 1].decision.tolist())
    neg_samples = np.array(score_df[score_df['status'] == 0].decision.tolist())

    mean_diff = pos_samples.mean() - neg_samples.mean()

    n_pos, n_neg = len(pos_samples), len(neg_samples)

    var_pos = pos_samples.var()
    var_neg = neg_samples.var()

    sd_pool = np.sqrt((n_pos * var_pos + n_neg * var_neg) / (n_pos + n_neg))
    cohen_d = mean_diff / sd_pool
    return cohen_d


def plot_decision_function(score_df, partition, output_file):
    """
    Plots the decision function for a given partition (either 'train' or
    'test') and saves a figure to file.

    Arguments:
    :param score_df: a specific folds decision scores and status
    :param partition: either 'train' or 'test' will plot performance
    :param output_file: file to output the figure
    """
    ax = sns.kdeplot(score_df.ix[(score_df.status == 1) &
                                 (score_df.partition == partition), :]
                     .decision, color='red', label='Deficient',
                     shade=True)
    ax = sns.kdeplot(score_df.ix[(score_df.status == 0) &
                                 (score_df.partition == partition), :]
                     .decision, color='blue', label='Wild-Type',
                     shade=True)
    ax.set(xlabel='Decision Function', ylabel='Density')
    ax.set_title('Classifier Decision Function')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def interpolate_roc(df, by=0.01):
    """
    Interpolates true positive rates for a given set of false positive points.

    Usage: pandas.DataFrame.apply(interpolate_roc) on roc output df

    Arguments:
    :param df: pandas DataFrame with fpr and tpr data
    :param by: how many points to arrange by in np.arange()
    """
    x_post = np.arange(0, 1, by)
    x_pre = [0.0] + list(df.fpr) + [1.0]
    y_pre = [0.0] + list(df.tpr) + [1.0]
    y_post = np.interp(x_post, xp=x_pre, fp=y_pre, left=0, right=1)
    return pd.DataFrame.from_items([('fpr', x_post), ('tpr', y_post)])


if __name__ == '__main__':
    # Load Command Arguments
    parser = argparse.ArgumentParser()

    # Required parameters
    parser.add_argument("-g", "--gene", dest="gene",
                        help="which gene to prepare cross validation")
    parser.add_argument("-t", "--tissue", dest="tissue",
                        help="which X matrix to load")
    parser.add_argument("-o", "--output_file", dest="out_file",
                        help="location to write the results to")
    parser.add_argument("-c", "--classifier", dest="classifier",
                        help="the type of classifier to use")
    parser.add_argument("-e", "--y_origin_file", dest="y_origin_file",
                        help="Where the y matrices are pulled from",
                        default='data/Y/')
    parser.add_argument("-y", "--hyperparameters", dest="hyperparameters",
                        help="hyperparameters to optimize in cross validation")

    # Optional parameters
    parser.add_argument("-x", "--alt_params", dest="alt_params",
                        help="additional set of hyperparameters to optimize in"
                             "cross validation", default=0.1)
    parser.add_argument("-f", "--num_folds", dest="num_folds", default=5,
                        help="the number of cross validation fold intervals")
    parser.add_argument("-z", "--holdout", dest="holdout", default=False,
                        help="boolean of extracting holdout set or not")
    parser.add_argument("-p", "--fraction_holdout", dest="fraction_holdout",
                        help="Fraction of samples to be heldout", default=0.1)
    parser.add_argument("-m", "--min_samples", dest="min_samples",
                        help="the minimum number of samples with a mutation to"
                        "consider the tissue", default=5)
    parser.add_argument("-i", "--iterations", dest="iterations",
                        help="number of random seed shuffles", default=50)
    parser.add_argument("-n", "--num_mad_genes", dest="num_mad_genes",
                        help="number of mad genes to subset", default=8000)
    parser.add_argument("-d", "--x_origin_file", dest="x_origin_file",
                        help="Where the x matrices are pulled from",
                        default='data/X/normalized/')
    parser.add_argument("-r", "--plot_roc", dest="roc",
                        help="decide if an ROC plot should be output",
                        action='store_true')
    parser.add_argument("-v", "--validation_sub", dest='validation_sub',
                        action='store_true')
    parser.add_argument("-w", "--coef_weights", dest='coefficients',
                        action='store_true')
    parser.add_argument("-s", "--effect_size", dest='effect_size',
                        help="determine effect sizes of cv folds",
                        action='store_true')
    args = parser.parse_args()

    # Load Constants
    # Required parameters
    gene = args.gene
    tissue = args.tissue.replace(',', '-')
    out_file = args.out_file
    classifier = args.classifier
    y_base = args.y_origin_file
    c = args.hyperparameters.split(',')

    # Optional parameters
    l = args.alt_params.split(',')
    num_folds = args.num_folds
    hold = args.holdout
    fraction_holdout = args.fraction_holdout
    min_samples = args.min_samples
    k = args.iterations
    random_seeds = random.sample(range(1, 1000000), int(k))
    num_mad_genes = args.num_mad_genes
    x_base = args.x_origin_file
    roc = args.roc
    validation_subset = args.validation_sub
    effect_size = args.effect_size

    # Generate filenames to save output plots
    fig_base = out_file.split('/')[1].split('.')[0]

    if roc:
        output_figure = os.path.join('figures',
                                     '{}_{}_{}alpha_{}_l1ratio_{}_roc.pdf'
                                     .format(tissue, classifier, fig_base,
                                             str(c[0]), str(l[0])))
        ensemble_roc = '{}_average_ensemble.tsv'.format(
            os.path.splitext(out_file)[0])

        if effect_size:
            cohens_d = []
            tot_dec_s = pd.DataFrame(columns=['decision', 'status',
                                              'partition'])
            decision_plot_train = os.path.join('figures',
                                               '{}_decision_plot_train.pdf'
                                               .format(fig_base))
            decision_plot_test = os.path.join('figures',
                                              '{}_decision_plot_test.pdf'
                                              .format(fig_base))
            cohen_plot = os.path.join('figures',
                                      '{}_cohens_d_plot.pdf'.format(fig_base))
            cohen_file = os.path.join('results',
                                      '{}_cohens_d.tsv'.format(fig_base))
    else:
        output_figure = os.path.join('figures',
                                     '{}_{}_{}_cv.png'.format(tissue,
                                                              classifier,
                                                              fig_base))
    # Are we interested in the weights/coefficients/genes?
    get_coef = args.coefficients
    weight_file = os.path.join('results',
                               '{}_{}_{}_weights.tsv'.format(tissue,
                                                             classifier,
                                                             fig_base))

    # Load and Process Data
    # RNAseq dictionary with sample IDs in order of columns
    sample_dictionary = pd.read_pickle(
                           os.path.join('output',
                                        'RNAseq_sample_dictionary.p'))

    # Generate Y Dictionary
    y_dict = generate_y_dict(sample_dictionary, gene, tissue, y_base)

    # Load X matrix
    if x_base == 'data/X/normalized/':
        x_fh = os.path.join(x_base,  '{}_X_normalized.tsv'.format(tissue))
    else:
        x_fh = os.path.join(x_base, '{}.tsv'.format(tissue))

    X = pd.read_csv(x_fh, delimiter='\t', index_col=0)

    # subset using MAD genes
    mad = pd.Series.from_csv(os.path.join('tables', 'full_mad_genes.tsv'),
                             sep='\t')
    mad = mad[:num_mad_genes]
    X = X.ix[mad.index]
    X = X.dropna()

    # Build a DataFrame of zeros for storing weights if option is selected
    empty_zero = np.zeros(X.shape[0])
    weight_coef = pd.DataFrame(empty_zero, index=X.index, columns=['weight'])

    if validation_subset:
        V = pd.read_csv(os.path.join('data', 'validation', 'normalized',
                                     'validation_set.tsv'),
                        delimiter='\t', index_col=0)
        V = V.ix[mad.index]
        V = V.dropna()
        X = X.ix[V.index]

    # roc output is iteratively made if ROC curve is to be output
    if roc:
        col_names = ['tpr', 'fpr', 'auc', 'type']
        roc_output = pd.DataFrame(columns=col_names)

    # Generate Cross Validation Parameters
    # Do this for each random seed
    cv_accuracy = []
    for seed in random_seeds:
        print(seed)
        # Build a cross validation partition object
        pancan_cv = cancer_cv(y_dict, num_folds, hold, min_samples, seed)

        # Determine fraction to holdout and assign folds for y and x matrices
        if hold:
            pancan_cv.holdout_samples(fraction_holdout)

        # Ensures equal partitioning into training and testing datasets
        pancan_cv.assign_folds()
        pancan_cv.assign_x_matrix_folds(X)

        # Perform num_folds cross validation on the total dataset
        for fold in range(0, num_folds):

            # Split according to the testing fold
            x_train, x_test, y_train, y_test = pancan_cv.split_train_test(fold)

            for param in c:
                if classifier == 'SVM':
                    clf = SVC(C=float(param), class_weight='balanced')
                    tr_err, te_err = train_test_error(clf, x_train, y_train,
                                                      x_test, y_test)
                    cv_accuracy.append(['train', tr_err, param, tissue, seed])
                    cv_accuracy.append(['test', te_err, param, tissue, seed])

                elif classifier == 'elasticnet':
                    for param_b in l:
                        clf = SGDClassifier(loss='log',
                                            penalty='elasticnet',
                                            alpha=float(param),
                                            l1_ratio=float(param_b),
                                            class_weight='balanced',
                                            random_state=seed)
                        if roc:
                            clf.fit(x_train, y_train)
                            y_score = clf.decision_function(x_test)
                            y_score_train = clf.decision_function(x_train)

                            fpr, tpr, _ = roc_curve(y_test, y_score,
                                                    drop_intermediate=False)
                            fpr_, tpr_, _ = roc_curve(y_train, y_score_train,
                                                      drop_intermediate=False)

                            test_auc = roc_auc_score(y_test, y_score,
                                                     average='weighted')
                            train_auc = roc_auc_score(y_train, y_score_train,
                                                      average='weighted')

                            temp_test = pd.DataFrame([tpr, fpr,
                                                     [test_auc] * len(tpr),
                                                     ['test'] * len(tpr)],
                                                     index=col_names).T
                            temp_train = pd.DataFrame([tpr_, fpr_,
                                                      [train_auc] * len(tpr_),
                                                      ['train'] * len(tpr_)],
                                                      index=col_names).T

                            full = temp_test.append(temp_train,
                                                    ignore_index=True)
                            full = full.assign(seed=seed)
                            full = full.assign(fold=fold)
                            roc_output = roc_output.append(full,
                                                           ignore_index=True)
                            if get_coef:
                                p_add = pd.DataFrame(clf.coef_.T,
                                                     index=x_test.columns,
                                                     columns=['weight'])
                                weight_coef = weight_coef.add(p_add,
                                                              fill_value=0)
                            if effect_size:
                                # Build input dataframes
                                y_st_df = pd.DataFrame(y_score_train,
                                                       index=y_train.index,
                                                       columns=['decision'])
                                y_st_df = y_st_df.assign(status=y_train)
                                y_st_df = y_st_df.assign(partition='train')

                                y_stes_df = pd.DataFrame(y_score,
                                                         index=y_test.index,
                                                         columns=['decision'])
                                y_stes_df = y_stes_df.assign(status=y_test)
                                y_stes_df = y_stes_df.assign(partition='test')

                                # Get effect size for training and testing
                                # separately, which will quantify the
                                # magnitude of separation between groups
                                train_cohen = get_cohens_d(y_st_df)
                                test_cohen = get_cohens_d(y_stes_df)

                                # Compile results
                                cohens_d.append(['train', train_cohen])
                                cohens_d.append(['test', test_cohen])
                                tot_dec_s = tot_dec_s.append(y_stes_df)
                                tot_dec_s = tot_dec_s.append(y_st_df)

                        else:

                            train_err, test_err = train_test_error(clf,
                                                                   x_train,
                                                                   y_train,
                                                                   x_test,
                                                                   y_test)
                            cv_accuracy.append(['train', train_err, param,
                                                param_b, tissue, seed])
                            cv_accuracy.append(['test', test_err, param,
                                                param_b, tissue, seed])

    # Write results to file
    if roc:
        roc_output.to_csv(out_file, sep='\t', index=False)

        # Interpolate ROC curve and aggregate for average ROC
        interpolated_df = roc_output.groupby(['fold', 'seed', 'type']) \
                                    .apply(interpolate_roc) \
                                    .reset_index() \
                                    .drop('level_3', 1)
        avg_df = interpolated_df.groupby(['type', 'fpr']).tpr.mean() \
                                                             .reset_index()

        zero_df = pd.DataFrame.from_items([('type', ['test', 'train'] * 2),
                                           ('fpr', [0, 0, 1, 1]),
                                           ('tpr', [0, 0, 1, 1])])
        avg_df = avg_df.append(zero_df, ignore_index=True)
        avg_df = avg_df.sort_values(['type', 'tpr']).reset_index(drop=True)
        avg_df.to_csv(ensemble_roc, sep='\t', index=False)
    else:
        cv_accuracy = pd.DataFrame(cv_accuracy)
        if classifier == 'SVM':
            cv_accuracy.columns = ['class', 'error', 'param', 'tissue', 'seed']
        elif classifier == 'elasticnet':
            cv_accuracy.columns = ['class', 'error', 'alpha', 'l1_ratio',
                                   'tissue', 'seed']
        cv_accuracy.to_csv(out_file, index=False, sep='\t')

    if get_coef:
        weight_coef['abs'] = weight_coef['weight'].abs()
        sorted_weight = weight_coef.sort_values('abs', ascending=False)
        sorted_weight = sorted_weight['weight']
        sorted_weight.to_csv(weight_file, index=True, sep='\t')

        command = 'R --no-save --args ' + weight_file + ' ' + \
                  'figures/weight_genes.pdf < scripts/viz/viz_geneweights.r'
        call(command, shell=True)

    # Plot
    if roc:
        command = 'R --no-save --args ' + out_file + ' ' + ensemble_roc + \
                  ' ' + output_figure + ' < ' + 'scripts/viz/viz_roc.r'
    else:
        command = 'R --no-save --args ' + out_file + ' ' + output_figure + \
                  ' ' + classifier + ' < ' + 'scripts/viz/viz_cv.r'
    call(command, shell=True)

    if effect_size:
        sns.set(font_scale=2)
        sns.set_style("white")

        # Plot Decision Functions
        plot_decision_function(tot_dec_s, 'train', decision_plot_train)
        plot_decision_function(tot_dec_s, 'test', decision_plot_test)

        # Violin Plots of Cohen's D Effect sizes
        cohens_df = pd.DataFrame(cohens_d, columns=['partition', 'cohens_D'])
        ax = sns.violinplot(x='partition', y='cohens_D', data=cohens_df,
                            palette="Set2")
        ax.set(xlabel='', ylabel="Cohen's D")
        ax.set_title('Cross Validation Effect Sizes')
        sns.despine()
        plt.tight_layout()
        plt.savefig(cohen_plot)
        plt.close()

        # Save Cohen's D effect size estimates
        cohens_df.to_csv(cohen_file, sep='\t', index=False)
