FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-mosdepth/bin/:/root/google-cloud-sdk/bin/:${PATH}

RUN apt-get -y update && \
    apt-get -y install curl zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev && \
    apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# activate conda environment
RUN echo "source activate lr-mosdepth" > ~/.bashrc

# copy the script
COPY coverage_stats.py /
