FROM ubuntu:20.04
LABEL maintainer="https://github.com/twolinin/longphase"
LABEL version="1.5.1"

RUN apt-get update && \
    apt-get install -y git g++ gcc autoconf make zlib1g-dev libbz2-dev liblzma-dev && \
    rm -rf /var/lib/apt/lists/* 

WORKDIR /opt/longphase
RUN git clone https://github.com/twolinin/longphase.git /opt/longphase && \
    autoreconf -i && \
    ./configure && \
    make -j 4 && \
    rm -rf /opt/longphase/.git

ENV PATH="${PATH}":${HOME}/bin:/opt/longphase
