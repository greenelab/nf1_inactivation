"""
(C) Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/build_y.py

Determine y matrix for each tissue according to mutation status and gene id

Usage:
Run on the command line 'build_tissue_specific_y_matrices.py'

    With required flags:

    --gene      The gene to determine mutation status for (e.g. 'NF1')
    --tissue    Determine filter for specific subset of tissues (default 'GBM')

    And the following optional flags:

    --filter     Filter mutation types from gold standard (default 'Silent')
    --unclass    RNAseq samples without mutations are filtered
    --y-origin   Location of where to save the Y matrices (default 'data/Y/')

Output:
Y matrix of mutation status for the input gene
"""

import pandas as pd
import argparse
from util.knn import get_y_class

####################################
# Load Command Arguments
####################################
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gene", dest="gene",
                    help="which gene to determine mut status")
parser.add_argument("-t", "--tissue", dest="tissue",
                    help="tissues to subset", default='GBM')
parser.add_argument("-f", "--filter", dest="filter",
                    help="classes of mutations to filter from gold standard",
                    default='Silent')
parser.add_argument("-u", "--unclassified", dest='unclass',
                    action='store_true', help="determine to remove samples "
                    "that have not been sequenced")
parser.add_argument("-y", "--y-origin", dest="out_fh",
                    help="Location of where to save the y files",
                    default='data/Y/')
args = parser.parse_args()

####################################
# Load Constants
####################################
GENE = args.gene
TISSUE = args.tissue
FILTER = args.filter
UNCLASSIFIED = args.unclass
OUT_FH = args.out_fh

####################################
# Load Data
####################################
mutation = pd.read_csv('data/PANCAN_mutation', delimiter='\t')

####################################
# Build Y matrices
####################################
# Generate the file name
file_name = OUT_FH + 'Y_' + TISSUE + '_' + GENE + '.tsv'

# Retrieve information about gold standard samples
mut_status, gold_neg, gold_pos, unclass, X = get_y_class(mutation, TISSUE,
                                                         [GENE], [FILTER])

# Generate dataframe
output_y = pd.DataFrame([mut_status, X.columns.values.tolist()]).T
output_y.columns = ['status', 'sample_id']

# If we decide to remove unsequenced samples, build Y w/o
if UNCLASSIFIED:
    output_y = output_y[~output_y['sample_id'].isin(unclass)]

# Generate file name and output to file
output_y.to_csv(file_name, sep='\t', index=False)
