# Modified Dockerfile script from other long-read-pipelines scripts to download bowtie2 and samtools

# Start with a good base python3 image:
FROM ubuntu:20.04
MAINTAINER Raphael Brosula <brosula@broadinstitute.org>

# Make sure we don't need to interact with any package installations:
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /
WORKDIR /

########################################################################################################################
# DEPENDENCIES
RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y build-essential

# Dependencies for samtools:
RUN apt-get install -y bzip2 curl gnupg2 libc-dev ncurses-dev libcurl4-openssl-dev libssl-dev libbz2-dev liblzma-dev zlib1g-dev

# Additional Dependencies:
RUN apt install -y curl wget pkg-config zip unzip default-jre
########################################################################################################################
# SOFTWARE:

RUN mkdir -p /usr/local/bin /usr/local/lib /usr/local/etc

# Samtools:
# Get samtools source:
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar -xjf samtools-1.11.tar.bz2 && \
	cd samtools-1.11 && \
	./configure && make install && \
    cd .. && \
    rm -rf samtools-1.11 samtools-1.11.tar.bz2

# bowtie2/2.5.3
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.3/bowtie2-2.5.3-linux-x86_64.zip && \
    unzip bowtie2-2.5.3-linux-x86_64.zip && \
    cd bowtie2-2.5.3-linux-x86_64 && \
    mv bowtie2* /usr/local/bin && \
    mkdir -p /usr/local/etc/bowtie2 && \
    mv AUTHORS BOWTIE2_VERSION LICENSE MANUAL MANUAL.markdown NEWS README.md /usr/local/etc/bowtie2 && \
    cd .. \
    rm -rf bowtie2-2.5.3-linux-x86_64.zip bowtie2-2.5.3-linux-x86_64/

########################################################################################################################
########################################################################################################################
########################################################################################################################

# Other utilities:
RUN apt-get clean
RUN apt install -y vim emacs nano
RUN apt-get clean