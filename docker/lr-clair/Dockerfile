FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/clair/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y git wget bzip2 make gcc vcftools && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/clair

# install anaconda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh && \
    bash Anaconda3-2019.10-Linux-x86_64.sh -b -p /opt/conda && \
    rm Anaconda3-2019.10-Linux-x86_64.sh

# create conda environment
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clair-env -c bioconda -y clair
RUN echo "source activate clair-env" > ~/.bashrc
ENV PATH /opt/conda/envs/clair-env/bin:$PATH
RUN /bin/bash -c ". activate clair-env && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install --no-cache-dir intervaltree"

# download the trained model for ONT
RUN mkdir -p models/ont && \
    cd models/ont && \
    wget http://www.bio8.cs.hku.hk/clair_models/ont/122HD34.tar && \
    tar -xf 122HD34.tar

# download the trained model for PacBio CCS
RUN mkdir -p models/pacbio && \
    cd models/pacbio && \
    wget http://www.bio8.cs.hku.hk/clair_models/pacbio/ccs/15.tar && \
    tar -xf 15.tar

# download the trained model for Illumina
RUN mkdir -p models/illumina && \
    cd models/illumina && \
    wget http://www.bio8.cs.hku.hk/clair_models/illumina/12345.tar && \
    tar -xf 12345.tar

# clone Clair
RUN git clone --depth 1 https://github.com/HKU-BAL/Clair.git
ENV CLAIR=/opt/conda/envs/clair-env/bin/clair.py
