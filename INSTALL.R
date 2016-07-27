# (C) Gregory Way 2016
# NF1 Inactivation Classifier for Glioblastoma
# INSTALL.R
#
# Usage:
# Run by ANALYSIS.sh to reproduce compute environment

library('methods')

mirror <- "http://cran.us.r-project.org"
install.packages("checkpoint", repos = mirror)

library("checkpoint")

#####################
# Install CRAN packages
#####################
checkpoint("2016-07-25")

######################
# Install bioconductor packages
######################
source("https://bioconductor.org/biocLite.R")
bioc_pkgs <- c(
  'pd.hta.2.0',
  'hta20stprobeset.db',
  'hta20sttranscriptcluster.db',
  'annotate',
  'affxparser',
  'oligo',
  'WGCNA'
)

biocLite(bioc_pkgs, suppressUpdates = TRUE)

######################
# Install source packages
######################
# https://github.com/greenelab/TDM
devtools::install_github("greenelab/TDM")

