"""
Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/validation.py

Applies classifier with optimal hyperparameters to a validation dataset

Usage:
The script is run by run_pipeline.sh will apply classifier to a validation set

    With required flags:

    --validation_file    The file location of the validation set
    --x-matrix           The file location for the training data
    --out_file           The location to save the predictions
    --hyperparameters    The optimal alpha value
    --alt_params         The optimal l1 ratio for elastic net

    And optional flags:

    --ensemble           Number of iterations over training and testing folds
    --gene               Gene to create y dictionaries from (defaults to 'NF1')
    --tissue             Tissue to extract (defaults to 'GBM')
    --num_mad_genes      Number of MAD genes to subset (defaults to 8000)
    --y_origin_file      Where y matrices are stored (defaults to 'data/Y/')
    --num_folds          The number of cross validation folds (defaults to 5)
    --min_samples        Mimimum number of samples to consider (defaults to 5)
"""

import os
import argparse
import pandas as pd
import random
from subprocess import call
from process_rnaseq import normalize_data
from nf1_classifier import generate_y_dict
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import roc_auc_score
from util.cancer_cv import cancer_cv

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--validation_file", dest="validation_file",
                    help="file handle for validation data")
parser.add_argument("-x", "--x_matrix", dest="training_data",
                    help="file handle for training data")
parser.add_argument("-o", "--out_file", dest="out_file",
                    help="file handle for predictions")
parser.add_argument("-a", "--hyperparameters", dest="hyperparameters",
                    help="hyperparameters to optimize using cross validation")
parser.add_argument("-l", "--alt_params", dest="alt_hyperparameters",
                    help="additional set of hyperparameters to optimize using"
                         "cross validation", default=0.1)

parser.add_argument("-e", "--ensemble", dest="num_seeds",
                    help="number of ensemble models", default=10)
parser.add_argument("-g", "--gene", dest="gene",
                    help="which gene to prepare cross validation",
                    default='NF1')
parser.add_argument("-t", "--tissue", dest="tissue",
                    help="which X matrix to load", default="GBM")
parser.add_argument("-n", "--num_mad_genes", dest="num_mad_genes",
                    help="number of mad genes to subset", default=8000)
parser.add_argument("-y", "--y_origin_file", dest="y_origin_file",
                    help="Where the y matrices are pulled from",
                    default='data/Y/')
parser.add_argument("-f", "--num_folds", dest="num_folds", default=5,
                    help="the number of cross validation fold intervals")
parser.add_argument("-m", "--min_samples", dest="min_samples",
                    help="the minimum number of samples with a mutation to"
                    "consider the tissue", default=5)
args = parser.parse_args()

# Load Constants
validation_file = args.validation_file
x_file = args.training_data
out_file = args.out_file
c = float(args.hyperparameters)
l = float(args.alt_hyperparameters)
k = args.num_seeds
random_seeds = random.sample(range(1, 1000000), int(k))
gene = args.gene
tissue = args.tissue.replace(',', '-')
num_mad_genes = args.num_mad_genes
num_folds = args.num_folds
y_base = args.y_origin_file
min_samples = args.min_samples
hold = False

# Generate filenames to save output plots
fig_file = os.path.join('figures', out_file.split('/')[1].split('.')[0])
protein_data_file = os.path.join('data', 'protein_quantification.csv')

random.seed(123)

# Read Data
V = pd.read_csv(validation_file, delimiter='\t', index_col=0)
V = normalize_data(V, out_file=None, output=False, method='zscore')
X = pd.read_csv(x_file, delimiter='\t', index_col=0)

# subset using MAD genes
mad = pd.Series.from_csv(os.path.join('tables', 'full_mad_genes.tsv'),
                         sep='\t')
mad = mad[:num_mad_genes]
V = V.ix[mad.index].dropna()
X = X.ix[V.index]
samples = pd.DataFrame(V.columns.tolist()).T

# Analysis
if __name__ == '__main__':
    # RNAseq dictionary with sample IDs in order of columns
    sample_dictionary = pd.read_pickle(
                          os.path.join('output', 'RNAseq_sample_dictionary.p'))

    # Generate Y Dictionary
    y_dict = generate_y_dict(sample_dictionary, gene, tissue, y_base)

    col_names = ['sample_id', 'neg_prob', 'pos_prob', 'prediction', 'test_auc',
                 'seed']
    results_output = pd.DataFrame(columns=col_names)

    for seed in random_seeds:
        print(seed)
        seed_df = pd.DataFrame([seed] * samples.shape[1])

        # Build a cross validation partition object
        pancan_cv = cancer_cv(y_dict, num_folds, hold, min_samples, seed)

        # Ensures equal partitioning into training and testing datasets
        pancan_cv.assign_folds()
        pancan_cv.assign_x_matrix_folds(X)

        # Perform NUM_FOLDS cross validation on the total dataset
        for fold in range(0, num_folds):

            # Split according to the testing fold
            x_train, x_test, y_train, y_test = pancan_cv.split_train_test(fold)

            # Build classifier
            clf = SGDClassifier(loss='log', penalty='elasticnet', alpha=c,
                                l1_ratio=l, class_weight='balanced',
                                random_state=seed)

            # Fit the model using only the training set
            model = clf.fit(x_train, y_train)

            # Get the AUROC of test set
            y_score = model.decision_function(x_test)
            test_auc = roc_auc_score(y_test, y_score, average='weighted')
            test_auc = pd.DataFrame([test_auc] * samples.shape[1])

            # Also predict using only the training set
            predictions = pd.DataFrame(model.predict(V.T))
            probabilities = pd.DataFrame(model.predict_proba(V.T))
            results = pd.concat([samples.T, probabilities, predictions,
                                 test_auc, seed_df], axis=1)
            results.columns = col_names
            results_output = results_output.append(results, ignore_index=True)

    results_output.to_csv(out_file, index=False, sep='\t')

    # Visualize
    command = 'R --no-save --args ' + out_file + ' ' + protein_data_file + \
              ' ' + fig_file + ' < ' + 'scripts/viz/viz_validation.r'
    call(command, shell=True)
