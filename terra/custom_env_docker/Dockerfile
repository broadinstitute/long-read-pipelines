############### stage 0: build samtools and bcftools from source
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.18
ARG BCFTOOLS_VERSION=1.18
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
                 make \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}* && \
    apt-get -qqy purge autoconf automake bzip2 gcc make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help && \
    bcftools --help

############### stage 1: other commonly used bioinformatics utilities
FROM us.gcr.io/broad-dsde-methods/terra-jupyter-minimal-base:0.0.1

USER root

# copy from previous stage the binaries from samtools build
COPY --from=0 /usr/local/bin/* /usr/local/bin/

ARG lrmaCU="https://github.com/broadinstitute/lrma-cloud-utils.git"
ARG MINIMAP2_VERSION='2.26'
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                apt-transport-https \
                ca-certificates \
                gnupg \
                zlib1g-dev \
                aria2 \
                bc \
                bedtools \
                bzip2 \
                cron \
                curl \
                datamash \
                gawk \
                less \
                python3.10 \
                python3-pip \
                tabix \
                tmux \
                tree \
                vcftools \
                wget && \
    gcloud components update --quiet && \
    cd /tmp && \
    curl -L https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 | tar -jxvf - && \
    cd minimap2-${MINIMAP2_VERSION}_x64-linux && \
    mv k8 minimap2 paftools.js /usr/local/bin/ && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
        /var/tmp/* \
        /var/cache/apt/*

# key library
RUN python3 -m pip install --upgrade pip setuptools ipynbname && \
    python3 -m pip install --upgrade \
    git+${lrmaCU}

# terminal-customization
RUN cd /home/jupyter/ && \
    git clone \
        --depth=1 \
        https://github.com/Bash-it/bash-it.git /home/jupyter/.bash_it && \
    cd /home/jupyter/.bash_it && \
    bash install.sh --silent && \
    cd -

# cron-related
RUN echo 'jupyter' >> /etc/cron.allow && \
    echo 'jupyter' >> /etc/at.allow && \
    chmod u+s /usr/sbin/cron

# sendgrid for email notification
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                        postfix libsasl2-modules mailutils && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
        /var/tmp/* \
        /var/cache/apt/*
