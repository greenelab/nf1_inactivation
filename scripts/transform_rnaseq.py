"""
(C) Gregory Way 2016
NF1 Inactivation Classifier for Glioblastoma
scripts/transform_rnaseq.py

Will transform training RNAseq data using TDM to coerce distributions to appear
as microarray distributions

Usage:
The script is run by ANALYSIS.sh

    With required flags:

    --microarray-fh    File handle for the reference microarray distribution
    --rnaseq-fh        File handle for the target RNAseq data
    --out-fh           File handle to write the TDM transformation file
    --normalization    The normalization used to scale

The script will call an R script that outputs a training matched expression
matrix that is overwritten by zero-one normalization.
"""

import argparse
from subprocess import call
import pandas as pd
from process_rnaseq import zero_one_normalize

####################################
# Load Command Arguments
####################################
parser = argparse.ArgumentParser()  # Load command line options
parser.add_argument("-m", "--microarray-fh", dest="microarray_fh",
                    help="file handle for microarray data")
parser.add_argument("-r", "--rnaseq-fh", dest="rnaseq_fh",
                    help="file handle for RNAseq data")
parser.add_argument("-o", "--out-fh", dest="out_fh",
                    help="file handle for normalized output file")
parser.add_argument("-n", "--normalization", dest="normalization",
                    help="the normalization method")
args = parser.parse_args()

####################################
# Load Constants
####################################
MICRO_FH = args.microarray_fh
RNASEQ_FH = args.rnaseq_fh
OUT_FH = args.out_fh
NORMALIZATION = args.normalization

####################################
# Analysis
####################################
command = 'R --no-save --args ' + MICRO_FH + ' ' + RNASEQ_FH + \
          ' ' + OUT_FH + ' < ' + 'scripts/util/transform.r'
call(command, shell=True)

X = pd.read_csv(OUT_FH, delimiter='\t', index_col=0)
zero_one_normalize(X, OUT_FH, method=NORMALIZATION)
