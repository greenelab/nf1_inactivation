"""
Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/transform_rnaseq.py

Will transform training RNAseq data using TDM to coerce distributions to appear
as microarray distributions

Usage:
The script is run by run_pipeline.sh

    With required flags:

    --microarray_file    File location for the reference microarray data
    --rnaseq_file        File location for the target RNAseq data
    --tdm_out_file       File location to write the TDM transformation file
    --normalization      The normalization used to scale

The script will call an R script that outputs a training matched expression
matrix that is overwritten by the given input normalization scheme.
"""

import argparse
from subprocess import call
import pandas as pd
from process_rnaseq import normalize_data

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--microarray_file", dest="microarray_file",
                    help="file name for microarray data")
parser.add_argument("-r", "--rnaseq_file", dest="rnaseq_file",
                    help="file name for RNAseq data")
parser.add_argument("-o", "--tdm_out_file", dest="tdm_out_file",
                    help="file name for normalized output file")
parser.add_argument("-n", "--normalization", dest="normalization",
                    help="the normalization method")
args = parser.parse_args()

# Load Constants
micro_file = args.microarray_file
rnaseq_file = args.rnaseq_file
tdm_out_file = args.tdm_out_file
normalization = args.normalization

# Analysis
command = 'R --no-save --args {} {} {} < \
          scripts/util/transform.R'.format(micro_file, rnaseq_file,
                                           tdm_out_file)
call(command, shell=True)

x_matrix_tdm = pd.read_csv(tdm_out_file, delimiter='\t', index_col=0)
normalize_data(x_matrix_tdm, tdm_out_file, method=normalization)
