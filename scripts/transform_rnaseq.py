"""
Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/transform_rnaseq.py

Will transform training RNAseq data using TDM to coerce distributions to appear
as microarray distributions

Usage:
The script is run by run_pipeline.sh

    With required flags:

    --microarray_file    File handle for the reference microarray distribution
    --rnaseq_file        File handle for the target RNAseq data
    --out_file           File handle to write the TDM transformation file
    --normalization      The normalization used to scale

The script will call an R script that outputs a training matched expression
matrix that is overwritten by zero-one normalization.
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
parser.add_argument("-o", "--out_file", dest="out_file",
                    help="file name for normalized output file")
parser.add_argument("-n", "--normalization", dest="normalization",
                    help="the normalization method")
args = parser.parse_args()

# Load Constants
micro_file = args.microarray_file
rnaseq_file = args.rnaseq_file
out_file = args.out_file
normalization = args.normalization

# Analysis
command = 'R --no-save --args {} {} {} < \
          scripts/util/transform.R'.format(micro_file, rnaseq_file, out_file)
call(command, shell=True)

X = pd.read_csv(out_file, delimiter='\t', index_col=0)
normalize_data(X, out_file, method=normalization)
