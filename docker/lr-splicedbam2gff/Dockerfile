FROM continuumio/miniconda3

MAINTAINER Jonn Smith <jonn@broadinstitute.org>

# Update apt:
RUN apt-get --allow-releaseinfo-change update 

# Install utilities:
RUN apt-get -y install wget vim tree

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

ENV PATH=/root/google-cloud-sdk/bin/:${PATH}

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# Install conda packages:
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install spliced_bam2gff
