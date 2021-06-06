FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

RUN apt-get -y update && \
    apt-get -y install curl zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev gcc make && \
    apt-get clean

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-whatshap/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# activate conda environment
RUN echo "source activate lr-whatshap" > ~/.bashrc
