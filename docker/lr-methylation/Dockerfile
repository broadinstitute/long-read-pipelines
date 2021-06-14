FROM nvidia/cuda:10.1-cudnn7-devel-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get -y install wget lsb-release libhdf5-dev zlib1g-dev lbzip2 gcc g++ make pciutils

ENV VERSION=v0.4
RUN wget "https://github.com/hasindu2008/f5c/releases/download/$VERSION/f5c-$VERSION-release.tar.gz" \
    && tar xvf f5c-$VERSION-release.tar.gz \
    && cd f5c-$VERSION/ \
    && ./scripts/install-hts.sh \
    && ./configure \
    && make cuda=1

ENV PATH=/f5c-$VERSION:$PATH
