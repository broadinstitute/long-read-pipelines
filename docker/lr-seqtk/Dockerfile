FROM ubuntu:20.04

ARG SEQTK_VERSION=1.3-1

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 ca-certificates \
                 datamash \
                 seqtk=${SEQTK_VERSION} \
                 zlib1g-dev && \
    apt-get -qqy purge wget make && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*


