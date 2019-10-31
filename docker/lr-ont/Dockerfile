FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-ont/bin:/miniconda/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install utilities
RUN apt-get -y update && apt-get install tree

# install gsutil
RUN curl https://sdk.cloud.google.com | bash
