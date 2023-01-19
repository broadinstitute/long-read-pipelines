FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-admixture/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-admixture/lib/

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3-pip liblzma-dev wget curl \
	&& apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install ADMIXTURE
RUN wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz \
	&& tar zxvf admixture_linux-1.3.0.tar.gz \
	&& mv dist/admixture_linux-1.3.0/admixture /usr/local/bin/

RUN echo "source activate lr-admixture" > ~/.bashrc
