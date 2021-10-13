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
RUN apt-get update 
RUN apt-get install -y apt-utils
RUN apt-get -y upgrade

# Development / prereqs for building software:
RUN apt-get -y install make gcc g++ autoconf sudo 

# Utilities / tools:
RUN apt-get -y install git bash vim time bc sed perl wget curl bzip2 man 

# Reqs for numpy/scipy:
RUN apt-get -y install libc-dev zlib1g-dev gfortran 

# Pysam remaining requirements:
RUN apt-get -y install liblzma-dev libbz2-dev

# Get libcurses:
RUN apt-get -y install libncurses5-dev libncursesw5-dev

# Get python 3.8:
RUN apt-get install -y software-properties-common
RUN apt-get install -y libssl-dev libncurses5-dev libsqlite3-dev libreadline-dev libgdm-dev libdb4o-cil-dev libpcap-dev libtk8.6
RUN cd /opt && \
		wget https://www.python.org/ftp/python/3.8.2/Python-3.8.2.tgz && \
		tar -xf Python-3.8.2.tgz && \
		cd Python-3.8.2 && \
		./configure --enable-optimizations && make altinstall && \
		python3.8 -m pip install --upgrade pip

# Install required Python packages: 
RUN python3.8 -m pip install cython
RUN python3.8 -m pip install mappy
RUN python3.8 -m pip install --no-cache-dir numpy
RUN python3.8 -m pip install --no-cache-dir matplotlib
RUN python3.8 -m pip install --no-cache-dir scipy
RUN python3.8 -m pip install pysam
RUN python3.8 -m pip install biopython

# Link python 3.8 to the global python executable:
RUN ln -s /usr/local/bin/python3.8 /usr/local/bin/python

################################################################################
# Install local settings:
RUN echo "set nowrap\nset number\nset hlsearch" > ~/.vimrc

################################################################################
# Install utils:

RUN python3.8 -m pip install pybedtools

RUN python3.8 -m pip install intervaltree kerneltree tqdm
RUN git clone https://github.com/BrooksLabUCSC/flair.git

RUN git clone https://github.com/lh3/minimap2.git \
    && cd minimap2 \
    && make \
    && mv minimap2 /usr/local/bin \
    && cd .. \
    && rm -rf minimap2

RUN apt-get update 
RUN apt-get install -y --fix-missing libcurl4-openssl-dev libssl-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
    && tar -xf samtools-1.10.tar.bz2 \
    && cd samtools-1.10 \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.10 \
    && rm -rf samtools-1.10.tar.bz2

# Add in requirements for FLAIR tools:
RUN apt-get install -y r-base 
RUN apt-get install -y r-base-dev r-base-core
RUN python3.8 -m pip install --no-cache-dir pandas 
RUN python3.8 -m pip install --no-cache-dir rpy2 

RUN apt-get install -y libxml2-dev

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R

# Get scanpy:
RUN python3.8 -m pip install scanpy

# Get SALMON:
ARG salmon_version=1.4.0
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${salmon_version}/salmon-${salmon_version}_linux_x86_64.tar.gz && \
    tar -zxf salmon-${salmon_version}_linux_x86_64.tar.gz &&  \
		mv salmon-latest_linux_x86_64 salmon-${salmon_version} &&  \
		mv salmon-${salmon_version}/bin/* /usr/local/bin/ && \
		mv salmon-${salmon_version}/lib/* /usr/local/lib/ && \
		rm -rf salmon-${salmon_version} &&  \
		rm salmon-${salmon_version}_linux_x86_64.tar.gz

# Get UMI-tools v1.1.1:
RUN pip install umi_tools==1.1.1

# Get the SSW package:
# NOTE: There is no versioning on this package...
RUN pip install ssw

# Clean up apt-get process:
RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install our scripts:
RUN mkdir /python_scripts 
COPY python/* /python_scripts/
RUN chmod +x /python_scripts/*.py

################################################################################
# Final runtime configuration:
# Let's start at the root:
WORKDIR /

