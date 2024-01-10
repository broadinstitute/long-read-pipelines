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

RUN apt-get update && apt-get -y upgrade

# Development / prereqs for building software:
# Utilities / tools:
# Get libcurses:
# install gsutil requirements:
RUN apt-get -y install make gcc g++ autoconf sudo && \
    apt-get -y install git bash vim time bc sed perl wget curl bzip2 man unzip && \
    apt-get -y install liblzma-dev libbz2-dev libncurses5-dev libncursesw5-dev && \
    apt-get --allow-releaseinfo-change update && \
    apt install -y curl git-lfs time datamash

# Install python:
RUN apt-get install -y python python3 python3-pip

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# install gsutil:
RUN curl https://sdk.cloud.google.com | bash

# Get Java:
# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-11-jdk && \
    apt-get install -y ant && \
    apt-get clean
    
# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
RUN export JAVA_HOME

# Get snPEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
	unzip snpEff_latest_core.zip

# Get tabix (for bgzip):
RUN apt-get install -y tabix

################################################################################
# Final runtime configuration:
# Let's start at the root:
WORKDIR /

