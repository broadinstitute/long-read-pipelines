FROM ubuntu:16.04

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3 python3-pip liblzma-dev wget curl \
	&& apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install translocator
RUN git clone https://github.com/HKU-BAL/Translocator.git \
    && cd Translocator/ \
    && git checkout 40bc3ea4b8178f7162bd9a0bbf839c0dc44891c9 \
    && git reset --hard \
    && mkdir build/ \
    && cd build/ \
    && cmake .. \
    && make

ENV PATH="/Translocator/bin/translocator-1.0.0/:/root/google-cloud-sdk/bin/:${PATH}"