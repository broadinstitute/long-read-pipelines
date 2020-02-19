FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "source activate lr-utils" > ~/.bashrc
ENV PATH=/opt/conda/envs/lr-utils/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt update
RUN apt install -y python2.7
RUN curl https://sdk.cloud.google.com | bash

# copy other resources
COPY detect_run_info.py /usr/local/bin/
COPY gatk.jar /usr/local/bin/
