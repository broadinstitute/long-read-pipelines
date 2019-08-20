# This image contains dependencies used by the cloud downloader system.

FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-cloud-downloader/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install python2.7
RUN apt update
RUN apt install -y python2.7

# install gsutil
RUN curl https://sdk.cloud.google.com | bash
