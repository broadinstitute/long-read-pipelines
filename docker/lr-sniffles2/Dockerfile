FROM continuumio/miniconda3

MAINTAINER Evie Wan

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-sniffles2/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-sniffles2/lib/

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3-pip liblzma-dev wget curl \
	&& apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash


RUN echo "source activate lr-sniffles2" > ~/.bashrc