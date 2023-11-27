FROM ubuntu:20.04
LABEL maintainer="https://github.com/twolinin/longphase"
LABEL version="1.5.1"

RUN apt-get update && \
    apt-get install -y git g++ gcc autoconf make zlib1g-dev libbz2-dev liblzma-dev wget && \
    rm -rf /var/lib/apt/lists/* 

WORKDIR /opt/longphase
RUN git clone https://github.com/twolinin/longphase.git /opt/longphase && \
    autoreconf -i && \
    ./configure && \
    make -j 4 && \
    rm -rf /opt/longphase/.git

ENV PATH="${PATH}":${HOME}/bin:/opt/longphase

ARG VERSION="v0.7.2"
ARG SAMTOOLS_VERSION=1.13
ARG BCFTOOLS_VERSION=1.13

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 g++ \
                 make \
                 wget \
                 zlib1g-dev

RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}*