FROM continuumio/miniconda3:22.11.1

MAINTAINER Jonn Smith

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# copy python scripts
COPY python/* /python_scripts/
RUN chmod +x /python_scripts/*

# activate conda environment
RUN echo "source activate lr-malaria" > ~/.bashrc
