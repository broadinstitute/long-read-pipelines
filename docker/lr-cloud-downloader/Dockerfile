FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-cloud-downloader/bin/:/root/google-cloud-sdk/bin/:${PATH}
RUN echo "source activate lr-cloud-downloader" > ~/.bashrc

# install gsutil
RUN apt update
RUN apt install -y curl git-lfs
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc aria2 python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

