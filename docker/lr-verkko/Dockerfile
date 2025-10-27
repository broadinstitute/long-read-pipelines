############### stage 0
FROM continuumio/miniconda3:25.3.1-1 AS build
ARG PKM='conda'
ARG CONDA_ENV_NAME='verkko'
ARG INSTALL_CMD="${PKM} env create -n ${CONDA_ENV_NAME}"
ARG CLEAN_CMD="${PKM} clean --all --yes"

# conda-pack is needed to create a standalone conda environment that can be moved between systems.
# It allows you to package all the dependencies of this project into a single tarball.
# This tarball can then be unpacked on any system that has conda installed,
# allowing you to quickly and easily recreate your environment without transfering over
# the large conda cache, which can be several GB in size. This is especially useful for
# saving space on a docker image.
# Install conda-pack:
RUN ${PKM} update -y -n base conda && \
    ${PKM} install -y -c conda-forge conda-pack libmamba && \
    ${PKM} config --set solver libmamba && \
    ${CLEAN_CMD}

# Copy environment.yml to the container:
COPY environment.yml .
RUN ${INSTALL_CMD}  \
      --file environment.yml && \
    ${CLEAN_CMD}

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack \
    -n ${CONDA_ENV_NAME} \
    -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack

# We now have to fix the verkko executable because they don't know the difference between 
# Bourne shell and BASH shell:
RUN sed -i '1c #!/usr/bin/env bash' /venv/bin/verkko


RUN ${CLEAN_CMD}

############### stage 1 
# Build Rukki:
FROM rust:1.90.0-trixie AS build_rukki

RUN git clone --recursive https://github.com/marbl/rukki.git && \
      cd rukki && git checkout 719569528805be1a5ab03c4870814c2e88e1100f && \
      cargo build --release && cargo install --path .

############### stage 2
# Build Meryl / Merqury
FROM ubuntu:22.04 AS build_merqury
ARG MERYL_VERSION=1.4.1
ARG SAMTOOLS_VERSION=1.22.1
ARG BCFTOOLS_VERSION=1.22

# Meryl:
RUN apt-get update && apt-get install -yqq build-essential git wget && \
    wget https://github.com/marbl/meryl/releases/download/v${MERYL_VERSION}/meryl-${MERYL_VERSION}.Linux-amd64.tar.xz && \
    tar -xJf meryl-${MERYL_VERSION}.Linux-amd64.tar.xz

# Merqury:
RUN git clone https://github.com/marbl/merqury.git && \
      cd merqury && git checkout 1ad7c32 && export MERQURY=$PWD

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 build-essential \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 libncursesw5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 g++ \
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
    rm -rf bcftools-${BCFTOOLS_VERSION}*

############### stage 3
FROM ubuntu:22.04 AS runtime

MAINTAINER Jonn Smith <jonn@broadinstitute.org>
LABEL org.opencontainers.image.authors="Jonn Smith <jonn@broadinstitute.org>"

#### Basic utilities
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 build-essential \
                 ca-certificates \
                 curl \
                 gnupg \
                 pigz \
                 vim \
                 pv \
                 tree \
                 zlib1g-dev \
                 openjdk-11-jre \
                 curl \
                 bc \
                 git-lfs && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* 

# Copy /venv from the previous stages:
COPY --from=build /venv /venv
COPY --from=build_rukki /usr/local/cargo/bin/rukki /usr/local/cargo/bin/rukki
COPY --from=build_merqury /meryl-*/bin/* /usr/local/bin/
COPY --from=build_merqury /meryl-*/lib/* /usr/local/lib/
COPY --from=build_merqury /merqury /merqury

COPY --from=build_merqury /usr/local/bin/* /usr/local/bin/

ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:/usr/local/cargo/bin:/root/google-cloud-sdk/bin:$PATH"
ENV MERQURY=/merqury

# Now install gsutil so we have it in this conda env later:
RUN export PATH="$VIRTUAL_ENV/bin:/usr/local/cargo/bin:$PATH" && curl https://sdk.cloud.google.com | bash

## Setup crcmodc for gsutil:
#RUN apt-get install -y gcc aria2 python3-dev python3-setuptools && \
#			pip3 uninstall -y crcmod && \
#			pip3 install --no-cache-dir -U crcmod
#

