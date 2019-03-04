# This base image contains dependencies used by the PacBio EAP processing pipeline.
# It uses bitnami/minideb:stretch as its base image to work with cromwell.
# Most of this is lifted from Ben Weisburd's STR analysis repository (https://github.com/bw2/str-callers).

FROM bitnami/minideb:stretch

MAINTAINER Kiran Garimella

ENV PYTHON_VERSION="3.7.1"
ENV SAMTOOLS_VERSION="1.9"

RUN install_packages ca-certificates wget pkg-config gcc g++ make bzip2 zlib1g-dev dpkg-dev build-essential libcurl4-openssl-dev libbz2-dev liblzma-dev libcrypto++-dev libncursesw5-dev libgdbm-dev libc6-dev zlib1g-dev libsqlite3-dev libssl-dev openssl libffi-dev vim less curl xterm

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure --without-curses \
	&& make \
	&& make install

# python2 utils
#RUN install_packages python-dev python-setuptools python-pip \
#	&& pip install --upgrade wheel

# install python3.7.1
#RUN apt-get update \
#	&& apt-get dist-upgrade -y \
#	&& apt-get install -y python-smbus libncursesw5-dev libgdbm-dev libc6-dev zlib1g-dev libsqlite3-dev libssl-dev openssl libffi-dev
#
#RUN wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
#	&& tar xvf Python-${PYTHON_VERSION}.tar.xz \
#	&& rm Python-${PYTHON_VERSION}.tar.xz \
#	&& cd Python-${PYTHON_VERSION} \
#	&& ./configure --enable-optimizations \
#	&& make install

# python3 packages
#RUN python3 -m pip install --upgrade pip \
#	&& python3 -m pip install --upgrade tqdm pyvcf pyfaidx cyvcf2

# miniconda installation
#RUN curl -LO https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
#RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
#RUN rm Miniconda3-latest-Linux-x86_64.sh
#ENV PATH=/root/google-cloud-sdk/bin/:/miniconda/bin:${PATH}
#RUN conda update -y conda

# configure conda
#RUN conda config --add channels defaults
#RUN conda config --add channels bioconda
#RUN conda config --add channels conda-forge
#RUN conda update --all

# install conda packages
#RUN conda install -c bioconda pbccs
#RUN conda install -c bioconda pbmm2
#RUN conda install -c bioconda pbsv
#RUN conda install -c bioconda nanoplot
#RUN conda install -c bioconda fastqc

# install gsutil
#RUN curl https://sdk.cloud.google.com | bash

# copy over other resources
COPY bashrc /root/.bashrc
