FROM ubuntu:22.04

MAINTAINER Jonn Smith

# Install required packages:
RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y --no-install-recommends \
	wget \
	make \
	g++ \
	gcc \
	git \
	vim \
	unzip \
	bzip2 \
	curl \
	gnupg2 \
	git-lfs \
	parallel \
	libc-dev \
	ncurses-dev \
	libcurl4-openssl-dev \
	libssl-dev \
	libbz2-dev \
	liblzma-dev \
	gcc \
	python3-dev \
	python3-setuptools \
	python3-pip 

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN pip3 uninstall -y crcmod && \
    pip3 install --no-cache-dir -U crcmod

# Install BWA-MEM2:
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    mv bwa-mem2-2.2.1_x64-linux /opt/ && \
    for f in $(find /opt/bwa-mem2-2.2.1_x64-linux/ -type f -name \*bwa-mem\* ) ; do ln -s $f /usr/local/bin/$(basename $f) ; done && \
    rm bwa-mem2-2.2.1_x64-linux.tar.bz2 

# Bowtie2:
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.3-linux-x86_64.zip && \
    mv bowtie2-2.3.4.3-linux-x86_64 /opt && \
    for f in $( find /opt/bowtie2-2.3.4.3-linux-x86_64/ -type f -name \*bowtie\* ) ; do ln -s $f /usr/local/bin/$(basename $f) ; done && \
    rm bowtie2-2.3.4.3-linux-x86_64.zip

# Install samtools 1.11:
#### Specific for google cloud support
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
		&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
		&& apt-get update -y \
		&& apt-get install google-cloud-sdk -y

# Get samtools source:
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
		&& tar -xjf samtools-1.11.tar.bz2 \
		&& cd samtools-1.11 \
		&& ./configure \
		&& make install
RUN rm -rf /samtools-1.11 /samtools-1.11.tar.bz2

# Bedtools:
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make install && \
    cd .. && \
    rm -rf bedtools2 bedtools-2.31.1.tar.gz

# Cleanup
RUN rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* \
           /root/.cache/pip \
           /root/.cache/pip3 \
           /var/lib/apt/lists/* && \
	apt-get -qqy clean && \
    apt-get -qqy autoclean && \
    apt-get -qqy autoremove

