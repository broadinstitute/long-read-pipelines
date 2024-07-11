FROM rocker/r-ubuntu:22.04

MAINTAINER Jonn Smith

########################################################################################################################

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# Update our packages:
RUN apt-get install -y --no-install-recommends \
      build-essential autoconf autoconf-archive libcppunit-dev
