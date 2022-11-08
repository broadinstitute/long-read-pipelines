# Start with a good base python3 image:
FROM ubuntu:20.04
MAINTAINER Jonn Smith <jonn@broadinstitute.org>

# Make sure we don't need to interact with any package installations:
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /
WORKDIR /

################################################################################
# Install system packages:

# Make sure we can actually talk to package repos:
RUN apt-get update && apt-get install -y apt-utils software-properties-common && apt-get -y upgrade && add-apt-repository ppa:deadsnakes/ppa && apt-get -y update
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake \
			gcc g++ zlib1g zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev vim wget unzip curl vim python3.7 python3.7-distutils \
			libc-dev ncurses-dev bzip2 gnupg2 libcurl4-openssl-dev libssl-dev bc \
          && rm -rf /var/lib/apt/lists/*

# Install google cloud support:
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
		&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
		&& apt-get update -y \
		&& apt-get install google-cloud-sdk -y

# Install samtools:
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
		&& tar -xjf samtools-1.11.tar.bz2 \
		&& cd samtools-1.11 \
		&& ./configure \
		&& make install \
		&& cd ../ && rm -rf /samtools-1.11 /samtools-1.11.tar.bz2

# Install python and some requirements:
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3.7 get-pip.py && \
    pip3 install pysam tqdm editdistance numpy

# Install ssw lib:
ADD ssw/ssw_lib.py /root/.local/lib/python3.7/site-packages/ssw/ssw_lib.py
RUN touch /root/.local/lib/python3.7/site-packages/ssw/__init__.py \
    && wget -O ssw-zip.zip https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/master.zip \
	&& unzip ssw-zip.zip \
	&& cd Complete-Striped-Smith-Waterman-Library-master/src \
	&& make libssw.so \
	&& cd ../.. \
    && mkdir -p /usr/local/lib/ssw \
	&& mv Complete-Striped-Smith-Waterman-Library-master/src/libssw.so /usr/local/lib/ssw/libssw.so \
	&& rm -rf Complete-Striped-Smith-Waterman-Library-master

# Install bamtools:
RUN git clone https://github.com/pezmaster31/bamtools.git \
    && cd bamtools \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && make install \
    && cd ../../ \
    && rm -rf bamtools

# Install our scripts:
RUN mkdir /scripts
COPY scripts/* /scripts/
RUN chmod +x /scripts/*.py /scripts/*.sh

################################################################################
# Final runtime configuration:
# Let's start at the root:
WORKDIR /

