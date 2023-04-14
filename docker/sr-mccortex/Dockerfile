FROM ubuntu:latest

RUN apt-get update \
    && apt-get install -y git autoconf automake pkg-config make gcc g++ liblzma-dev libbz2-dev libz-dev libcurl4-openssl-dev libncurses-dev \
    && mkdir -p /opt && git clone --recursive https://github.com/mcveanlab/mccortex /opt/mccortex \
    && rm -rf /var/lib/{apt,dpkg,cache,log}

WORKDIR /opt/mccortex
RUN sed -i -e 's/then autoreconf/then autoreconf --install/g' -e 's/--disable-lzma/--disable-lzma --host x86_64/' /opt/mccortex/libs/Makefile \
    && make MAXK=31 all && install -p  bin/* /usr/local/bin
