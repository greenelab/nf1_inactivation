# (C) Gregory Way 2016
# NF1 Inactivation Classifier for Glioblastoma
# INSTALL.R
#
# Usage:
# Run by the Docker build to reproduce compute environment

#####################
# Install CRAN packages
#####################
mirror <- "http://cran.us.r-project.org"

cran_pkgs <- c(
  'Rcpp',
  'readr',
  'codetools',
  'curl',
  'httr',
  'git2r',
  'plyr',
  'dplyr',
  'readr',
  'Rmisc',
  'calibrate',
  'data.table',
  'devtools'
)

install.packages(cran_pkgs, repos = mirror)

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

