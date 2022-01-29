FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-ref/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install packages
RUN apt-get update --allow-releaseinfo-change
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev curl git-lfs \
    && rm -rf /var/lib/apt/lists/*

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 \
    && tar xjf samtools-1.14.tar.bz2 \
    && rm samtools-1.14.tar.bz2 \
    && cd samtools-1.14 \
    && ./configure --without-curses \
    && make \
    && make install

# install TRF
RUN wget -O /usr/local/bin/trf https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 \
    && chmod 755 /usr/local/bin/trf
