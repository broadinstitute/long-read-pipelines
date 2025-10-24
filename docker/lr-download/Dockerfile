FROM continuumio/miniconda3

ARG DEBIAN_FRONTEND=noninteractive
ARG samtools_version=1.9

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 git \
                 ssh \
                 build-essential \
                 cmake \
                 vim \
                 apt-transport-https \
                 gnupg \
                 curl \
                 libssl-dev \
                 libffi-dev \
                 libreadline-dev \
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
                 zlib1g-dev \
                 unzip \
                 gzip

ENV PATH=/opt/conda/envs/lr-sv/bin/:/root/google-cloud-sdk/bin/:${PATH}
RUN curl https://sdk.cloud.google.com | bash 


RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    aws --version

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && make -s all all-htslib \
    && make install install-htslib \
    && rm -rf samtools-${samtools_version} \
    && samtools --help