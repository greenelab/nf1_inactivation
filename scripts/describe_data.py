'''
Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/describe_data.py

Get sample size and number of NF1 mutations

Usage:
Run on the command line 'python describe_data.py'

    There are also optional flags

    --type        Ignore different types of mutations (defaults to 'Silent')
    --tissue      Tissue to consider (defaults to 'GBM')
    --gene        The gene to query (defaults to 'NF1')
    --out-file    Output file to write mutation counts
    --mut-file    Mutation description table loc (tables/nf1_mutations.tsv)

Output:
Three distinct files including:
1) A tissue sample dictionary for RNAseq measured samples
2) A table tracking the number of mutations per tissue of specific genes
3) A table of the specific mutations for each of the NF1 mutated samples
'''

import pandas as pd
import pickle
import argparse

####################################
# Load Command Arguments
####################################
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--type", dest="mut_types",
                    help="different types of mutations to write",
                    default='Silent')
parser.add_argument("-t", "--tissue", dest="tissue",
                    help="tissues to query", default='GBM')
parser.add_argument("-g", "--gene", dest="gene",
                    help="genes to query", default='NF1')
parser.add_argument("-f", "--out-file", dest="out_fh",
                    help="filename of where to save output",
                    default='tables/sample_mut_table.tsv')
parser.add_argument("-o", "--mut-file", dest="mut_fh",
                    help="filename of where to save mutation table output",
                    default='tables/nf1_mutations.tsv')
args = parser.parse_args()

####################################
# Load Constants
####################################
MUT_TYPES = args.mut_types.split(',')
TISSUE = args.tissue
GENE = args.gene
OUT_FH = args.out_fh
MUT_FH = args.mut_fh

####################################
# Load data
####################################
mutation = pd.read_csv('data/PANCAN_mutation', delimiter='\t')
clinical = pd.read_csv('data/PANCAN_clinicalMatrix', delimiter='\t')
rnaseq = pd.read_csv('data/HiSeqV2', delimiter='\t', index_col=0)
tissue_dict = pd.read_csv('tables/tcga_dictionary.tsv', delimiter='\t')

####################################
# Determine Tissue Specific Samples
####################################
tissues = clinical._primary_disease.dropna().unique().tolist()
tissues.remove('FFPE pilot phase II')

# Get full tissue name
tissue_name = tissue_dict[tissue_dict['acronym'] == TISSUE]['tissue'].values[0]

# Initialize a tissue sample dictionary be used to store sample IDs
tissue_dict_samples = {}

# Subset clinical information to tissue
clin_sub = clinical.loc[clinical['_primary_disease'] == tissue_name, :]
samp_size = clin_sub.shape[0]

# What are the sample ids
sample_ids = clin_sub['sampleID'].tolist()

# RNAseq samples
rnaseq_samp = set(sample_ids) & set(rnaseq.columns)

# populate the tissue sample dictionary
tissue_dict_samples[TISSUE] = sample_ids

# Subset the mutation file to only these samples
mut_samples = mutation.loc[mutation['#sample'].isin(sample_ids), :]

# Subset the RNAseq file to only those with NF1 mutations
rnaseq_pos = set(mut_samples['#sample']) & set(rnaseq.columns)

# Get types of mutations
nf1_mutations = mut_samples.loc[mut_samples['gene'] == GENE, :]

# Filter mutations
for mut_type in MUT_TYPES:
    nf1_mutations = nf1_mutations.loc[mutation['effect'] != mut_type, :]

# Samples with non-silent mutations
nf1_mut_samples = set(nf1_mutations['#sample'].tolist())

# Number of samples with NF1 mutations and RNAseq measurements
nf1_rnaseq = nf1_mut_samples.intersection(rnaseq_samp)

# Compile info for input tissue
tissue_info = [
    TISSUE,
    tissue_name,
    samp_size,
    len(set(mut_samples['#sample'].tolist())),
    len(rnaseq_samp),
    len(rnaseq_pos),
    len(nf1_mut_samples),
    len(nf1_rnaseq),
]

# Subset the mutations we have RNAseq info for
rnaseq_nf1 = nf1_mutations[nf1_mutations['#sample'].isin(nf1_rnaseq)]

####################################
# Save output files to disk
####################################
with open('output/tissue_sample_dictionary.p', 'wb') as pickle_fh:
    pickle.dump(tissue_dict_samples, pickle_fh)

tissue_info = pd.DataFrame(tissue_info).T
tissue_info.columns = ['acronym', 'disease', 'total samples',
                       'samples with mutations', 'samples with RNAseq',
                       'samples with mutations and RNAseq',
                       'total NF1 mutations', 'NF1 mutations with RNAseq']
tissue_info.to_csv(OUT_FH, sep='\t', index=False)

use_cols = ['#sample', 'chr', 'start', 'end', 'reference', 'alt', 'gene',
            'effect', 'Amino_Acid_Change']
rnaseq_nf1_subset = rnaseq_nf1[use_cols]
rnaseq_nf1_subset.to_csv(MUT_FH, sep='\t', index=False)
