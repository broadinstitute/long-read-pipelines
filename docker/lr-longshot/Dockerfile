FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create scripts=3.6 -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-longshot/bin/:/root/google-cloud-sdk/bin/:${PATH}
RUN echo "source activate lr-longshot" > ~/.bashrc
RUN activate lr-longshot

# install packages
RUN apt update
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev curl git-lfs \
    && rm -rf /var/lib/apt/lists/*

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
    && tar xjf samtools-1.11.tar.bz2 \
    && rm samtools-1.11.tar.bz2 \
    && cd samtools-1.11 \
    && ./configure --without-curses \
    && make \
    && make install
