# Start with a good base python3 image:
FROM ubuntu:21.04
MAINTAINER Jonn Smith <jonn@broadinstitute.org>

# Make sure we don't need to interact with any package installations:
ARG DEBIAN_FRONTEND=noninteractive

## Install any needed packages specified in requirements.txt
#RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Set the working directory to /
WORKDIR /

################################################################################
# Install system packages:

# Make sure we can actually talk to package repos:
RUN apt-get update 

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

COPY tesserae-1.2.0.tar.gz /
RUN python3.8 -m pip install /tesserae-1.2.0.tar.gz
RUN rm /tesserae-1.2.0.tar.gz 

RUN git clone --recursive https://github.com/noporpoise/seq-align.git
RUN cd seq-align \
    && make

RUN git clone --recursive https://github.com/jonn-smith/MosaicAligner.git
RUN cd MosaicAligner \
    && make

RUN git clone --recursive https://github.com/jonn-smith/TALON.git
RUN cd TALON \
		&& python3.8 -m pip install .

RUN git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2 \
    && git checkout v2.29.2 \
    && make install

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
RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
    && tar -xf samtools-1.10.tar.bz2 \
    && cd samtools-1.10 \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.10 \
    && rm -rf samtools-1.10.tar.bz2

RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
		&& tar -xf bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
		&& rm bwa-mem2-2.0pre2_x64-linux.tar.bz2 

RUN git clone https://github.com/lh3/bwa.git \ 
		&& cd bwa \
		&& make

# Add in requirements for FLAIR tools:
RUN apt-get install -y r-base-dev r-base-core
RUN python3.8 -m pip install --no-cache-dir pandas 
RUN python3.8 -m pip install --no-cache-dir rpy2 

RUN apt-get install -y libxml2-dev

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R

# Install cartographer itself:
RUN mkdir /cartographer
COPY python/cartographer.py /cartographer/
COPY python/cartographer_bi_dir_segments.py /cartographer/
COPY python/extract_bounded_read_sections.py /cartographer/
COPY python/get_approx_raw_subread_array_lengths.py /cartographer/

COPY test_data/* /cartographer/

RUN chmod +x /cartographer/*.py
RUN ln -s /cartographer/cartographer.py /cartographer.py

COPY bash/quick_test_cartographer.sh /quick_test_cartographer.sh

################################################################################
# Copy the current directory contents into the container at /app
COPY pac_bio_left_adapter.fasta /reference_sequences/pac_bio_left_adapter.fasta

################################################################################
# Final runtime configuration:

# Clean up apt-get process:
RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Let's start at the root:
WORKDIR /

