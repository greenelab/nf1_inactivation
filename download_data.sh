#!/bin/bash                                                                                                                                                                           

##################
# NF1 Inactivation Classifier for Glioblastoma
# Gregory Way 2016
##################

# All data is publicly available and downloaded from UCSC Xena
# https://genome-cancer.ucsc.edu/proj/site/xena/datapages/?cohort=TCGA%20Pan-Cancer

# Because the database is continously updated and to ensure reproducibility,
# access data from zenodo cached download.

# RNAseq and Clincal data were downloaded on 8 March 2016
# Mutation data was downloaded on 12 June 2015

url=https://zenodo.org/record/56735/files/gbm_classifier_data.tar.gz
wget --directory-prefix data/ $url

# Extract and move files
tar -zxvf data/gbm_classifier_data.tar.gz
mv gbm_download/* data/
rm -rf gbm_download
rm data/gbm_classifier_data.tar.gz
