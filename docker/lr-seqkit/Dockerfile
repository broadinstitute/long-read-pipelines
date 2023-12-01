FROM ubuntu:20.04

ARG SEQKIT_VERSION=v2.4.0
ARG DL_LINK="https://github.com/shenwei356/seqkit/releases/download/${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz"

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 ca-certificates \
                 datamash \
                 wget \
                 zlib1g-dev && \
    cd /tmp && \
        wget ${DL_LINK} && tar -xf seqkit_linux_amd64.tar.gz && cp seqkit /usr/local/bin && \
    cd - && \
    apt-get -qqy purge wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

