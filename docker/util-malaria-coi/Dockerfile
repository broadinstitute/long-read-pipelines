FROM rocker/r-ubuntu:22.04

MAINTAINER Jonn Smith

########################################################################################################################

# Setup python3 and crcmodc for gsutil:
RUN apt-get update && \
    apt-get install -y gcc python3-dev python3-setuptools python3-pip && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

# Update our packages:
RUN apt-get install -y \ 
      build-essential autoconf autoconf-archive libcppunit-dev

# Install other dependencies / tools:
RUN apt-get install -y vim 

# Install DEploid:
RUN git clone https://github.com/DEploid-dev/DEploid.git && \ 
      cd DEploid && \
      ./bootstrap && \
      make install

# Install RealMcCOIL:
RUN git clone https://github.com/EPPIcenter/THEREALMcCOIL.git && \ 
    cd THEREALMcCOIL/categorical_method && \
    rm *.o *.so && \
    R CMD SHLIB McCOIL_categorical_code.c llfunction_het.c && \
    cd ../proportional_method && \
    rm *.o *.so && \
    R CMD SHLIB McCOIL_prop_code.c llfunction.c


