FROM continuumio/anaconda3:latest
MAINTAINER "Greg Way" <gregway@mail.med.upenn.edu>

# Install base packages
RUN apt-get update --fix-missing && apt-get install \
    libcurl4-openssl-dev \
    libssl-dev \
    libssh2-1-dev \
    vim

###########################
# Clone nf1_inactivation Repo
###########################
RUN hg clone https://gwaygenomics@bitbucket.org/gwaygenomics/nf1_inactivation

###########################
# Install R and Packages
###########################
RUN apt-get update --fix-missing && apt-get install \
    r-base \
    r-base-dev \
    r-cran-rcpp

RUN echo "deb http://cran.irsn.fr/bin/linux/debian jessie-cran3/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480
RUN apt-get update && apt-get install r-base

RUN Rscript nf1_inactivation/INSTALL.R

