FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /
COPY c3poa.config.txt /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-c3poa/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install git
RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev tree curl libncurses-dev \
	&& apt-get clean 

# install C3POa
RUN git clone https://github.com/rvolden/C3POa.git \
        && cd C3POa \
        && ./setup.sh

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# set LD library path
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-c3poa/lib/

# activate conda environment
RUN echo "source activate lr-c3poa" > ~/.bashrc
