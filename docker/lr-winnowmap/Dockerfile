FROM ubuntu:16.04

MAINTAINER Yulia Mostovoy

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3 python3-pip liblzma-dev wget curl \
	&& apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install translocator
RUN git clone https://github.com/marbl/Winnowmap.git \
    && cd Winnowmap/ \
    && git checkout 26a9ad1 \
    && git reset --hard \
    && make -j8

ENV PATH="/Winnowmap/bin:/root/google-cloud-sdk/bin/:${PATH}"