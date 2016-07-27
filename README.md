############################################
# NF1 Inactivation Classifier for Glioblastoma
#
## (C) 2016 Gregory Way, Robert Allaway, Rebecca Barnard, Stephanie Bouley, Brian O'Roak, Camilo Fadul, Yolanda Sanchez, and Casey Greene
############################################

###################
# SUMMARY
###################
The repository contains instructions to replicate and build upon a classifier trained to detect an NF1 inactivation signature in glioblastoma gene expression data. We leverage publicly available data from the Cancer Genome Atlas (TCGA) to tune a logistic classifier with an elastic net penalty using stochastic gradient descent. NF1 is a tumor suppressor that regulates RAS (a well characterized oncogene). When NF1 is inactivated, RAS signaling continues unabated leading to uncontrolled cell growth. Patients with neurofibromatosis type I (caused by heterozygous germline mutation of NF1) have a predisposition for multiple tumor types including optic gliomas, pheochromocytomas, and malignant peripheral nerve sheath tumors. Furthermore, NF1 is one of the most commonly mutated genes in glioblastoma. NF1 can be inactivated genetically or by other mechanisms including microRNAs or targeted degradation by the proteosome (McGillicuddy et al. 2009) . Therefore, detecting inactivation solely by sequencing the NF1 gene can result in false negatives. Because we have previously identified compounds that are synthetically lethal in NF1 inactivated cells (Wood et al. 2011), the ability to detect patients with NF1 inactivation signatures could inform treatment decisions.

###################
# CONTACT
###################
Please report all bugs and direct analysis questions to: `gregway@mail.med.upenn.edu`

Please direct all other correspondence to: `csgreene@mail.med.upenn.edu` or `yolanda.sanchez@dartmouth.edu`

###################
# ANALYSIS
###################
To reproduce and build upon our analyses, clone the repository and run `./ANALYSIS.sh`.

RECOMMENDED:
To ensure reproducibility our analyses should be performed using the provided [Docker image](https://hub.docker.com/r/gregway/nf1_inactivation).

###################
# DATA
###################
All data is publicly available. TCGA data used to train the classifier was retrieved from [UCSC Xena](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/). All of our validation data was deposited in [GEO](http://www.ncbi.nlm.nih.gov/geo/).

###################
# DEPENDENCIES
###################
All analyses are performed in Anaconda Python 3.5.1 (see `environment.yml`). Visualizations are built using R version 3.2.3. For specific package dependencies see `Docker/Dockerfile`

###################
# ACKNOWLEDGEMENTS
###################
This work was supported by the Genomics and Computational Biology graduate group at The Unversity of Pennsylvania (G.P.W); the Gordon and Betty Moore Foundation’s Data-Driven Discovery Initiative (grant number GBMF 4552 to C.S.G.); and the American Cancer Society (grant number IRG 8200327 to C.S.G.).
