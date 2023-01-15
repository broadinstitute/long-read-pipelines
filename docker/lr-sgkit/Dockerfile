FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "source activate lr-sgkit" > ~/.bashrc
ENV PATH=/opt/conda/envs/lr-sgkit/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git-lfs parallel
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod
