# Start with a good base python3 image:
FROM nvidia/cuda:10.1-devel-ubuntu18.04

MAINTAINER Jonn Smith <jonn@broadinstitute.org>

# Set the working directory to /
WORKDIR /

################################################################################
# Managed Package Installation:

# Make sure we can actually talk to package repos:
RUN apt-get update

# Development / prereqs for building software:
RUN apt-get -y install make gcc g++ python3 python3-pip autoconf sudo
#RUN pip3 install --upgrade pip

# Utilities / tools:
RUN apt-get -y install git bash vim time bc sed perl wget curl bzip2 man

# Reqs for numpy/scipy:
RUN apt-get -y install libc-dev zlib1g-dev python3-dev gfortran

# Pysam remaining requirements:
RUN apt-get -y install liblzma-dev libbz2-dev

# Clean up apt-get process:
RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install required Python packages:
RUN pip3 install mappy
RUN pip3 install --no-cache-dir numpy
RUN pip3 install --no-cache-dir matplotlib
RUN pip3 install --no-cache-dir scipy
RUN pip3 install pysam
RUN pip3 install tqdm

# Install a recent version of cmake:
WORKDIR /tmp
RUN mkdir /opt/cmake && \
    wget https://github.com/Kitware/CMake/releases/download/v3.14.4/cmake-3.14.4-Linux-x86_64.sh && \
    sh /tmp/cmake-3.14.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
RUN rm /tmp/cmake-3.14.4-Linux-x86_64.sh

################################################################################
# Manual Package Installation:

# Install shasta:
WORKDIR /opt
RUN git clone https://github.com/chanzuckerberg/shasta.git && \
    cd /opt/shasta && \
    git fetch && \
    git checkout $git_commit && \
    git submodule update --init && \
    rm -rf /opt/shasta/.git
WORKDIR /opt/shasta/build
RUN ../scripts/InstallPrerequisites-Ubuntu.sh
RUN cmake .. && make all && make install/strip && /opt/shasta/build/shasta-install/bin/shasta -h

# Install MarginPolish:
# Dependencies (many are already installed):
RUN apt-get -y install cmake make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev
RUN apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev
RUN apt-get -y install liblzma-dev libhdf5-dev
# MarginPolish:
RUN git clone https://github.com/UCSC-nanopore-cgl/marginPolish.git \
      && cd marginPolish \
      && git submodule update --init \
      && mkdir build \
      && cd build \
      && cmake .. \
      && make \
      && make install \
      && cp marginPolish /usr/local/bin/.

# HELEN Prereqs:
RUN apt-get -y install cmake python3 python3-dev
RUN pip3 install h5py tqdm numpy torchnet pyyaml

# HELEN needs PyTorch:
#RUN pip3 install https://download.pytorch.org/whl/cu100/torch-1.0.1.post2-cp36-cp36m-linux_x86_64.whl
RUN pip3 install torch torchvision
## HELEN:
RUN git clone https://github.com/kishwarshafin/helen.git \
      && cd helen \
      && git fetch \
      && git submodule update --init \
      && git pull origin master \
      && git checkout d7b4fed7858e1971b453f9743a7bda154acff323 \
      && ./build.sh

# SAMTools
RUN sudo apt-get -y install libncurses5-dev libncursesw5-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xjf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure && \
    make install

# JQ:
RUN apt-get install -y jq


################################################################################
# Copy the current directory contents into the container at /app
COPY pac_bio_left_adapter.fasta /reference_sequences/pac_bio_left_adapter.fasta

################################################################################
# Some useful tooling for vim:
RUN echo "set number\nset nowrap\nsyntax on\nset expandtab\n" > ~/.vimrc

################################################################################
# Final runtime configuration:

# Let's start at the root:
WORKDIR /


