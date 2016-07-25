'''
(C) Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/build_x.py

Load in normalized expression matrix and subset to tissue

Usage:
Is called by ANALYSIS.sh but can be run separately on the command line

        example: python build_x.py -t 'GBM'

        Where -t (--tissue) is comma separated list of tissue acronyms

Output:
1) Write tissue specific X matrices to data/ folder
2) Pickle dictionary of tissues and the order of the RNAseq samples in X
3) RNAseq sample size table in tables/
4) Descending order of MAD genes in tables/
'''

import pandas as pd
import pickle
import argparse

####################################
# Load Command Arguments
####################################
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tissue", dest="tissue",
                    help="tissues to write output files for")
args = parser.parse_args()

####################################
# Load Constants
####################################
TISSUE = args.tissue

###################
# Define Functions
###################


def write_expression_file(tissue_id, x_subset, loc='data/X/raw/'):
    '''
    writes the subset expression file to /data/ folder

    Arguments:
    :param tissue_id: the TCGA tissue acronym
    :param x_subset: a subset x matrix with tissue specific IDs
    :param loc: the parent folder to save file

    Return:
    nothing, writes to file
    '''
    file_name = loc + tissue_id + '_RNAseq_X_matrix.tsv'
    x_subset.to_csv(file_name, sep='\t')

###################
# Load data
###################
# Raw expression values
raw_express = pd.read_table('data/HiSeqV2', index_col=0)

# Tissue dictionary with sample IDs
sample_dictionary = pd.read_pickle('output/tissue_sample_dictionary.p')

###################
# Write subset X matrices to file
###################
rnaseq_dict = {}
sample_ids = sample_dictionary[TISSUE]
x_subset = raw_express.loc[:, raw_express.columns.isin(sample_ids)]
rnaseq_dict[TISSUE] = x_subset.columns.values.tolist()
write_expression_file(TISSUE, x_subset)

###################
# Save output file
###################
# RNAseq sample dictionary
with open('output/RNAseq_sample_dictionary.p', 'wb') as pickle_fh:
    pickle.dump(rnaseq_dict, pickle_fh)
