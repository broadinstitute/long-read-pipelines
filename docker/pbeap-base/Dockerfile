# This base image contains dependencies used by the PacBio EAP processing pipeline.
# It uses bitnami/minideb:stretch as its base image to work with cromwell.
# Most of this is lifted from Ben Weisburd's STR analysis repository (https://github.com/bw2/str-callers).

FROM bitnami/minideb:stretch

MAINTAINER Kiran Garimella

ENV PYTHON_VERSION="3.7.1"
ENV SAMTOOLS_VERSION="1.9"

RUN install_packages ca-certificates pkg-config gcc g++ make bzip2 zlib1g-dev dpkg-dev build-essential libcurl4-openssl-dev libbz2-dev liblzma-dev libcrypto++-dev libgdbm-dev libc6-dev zlib1g-dev libsqlite3-dev libssl-dev openssl libffi-dev
RUN install_packages python-pip python-dev
RUN install_packages wget curl vim tree time

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

# install gsutil
RUN curl https://sdk.cloud.google.com | bash
ENV PATH=/root/google-cloud-sdk/bin/:${PATH}

# copy over other resources
COPY bashrc /root/.bashrc
