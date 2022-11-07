FROM ubuntu:21.04
MAINTAINER Jonn Smith <jonn@broadinstitute.org>

# Make sure we don't need to interact with any package installations:
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /
WORKDIR /

################################################################################
# Install system packages:

# Make sure we can actually talk to package repos:
RUN apt-get update 

RUN apt-get install -y wget

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz && \
		tar -zxf ncbi-blast-2.13.0+-x64-linux.tar.gz && \
		rm ncbi-blast-2.13.0+-x64-linux.tar.gz

RUN apt-get install -y libgomp1 vim

ENV PATH="${PATH}:/ncbi-blast-2.13.0+/bin"

# Clean up apt-get process:
RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

