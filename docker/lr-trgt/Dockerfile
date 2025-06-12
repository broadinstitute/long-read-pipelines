FROM ubuntu:22.04

RUN apt-get update && \
        apt-get install -y wget bzip2 gcc g++ zlib1g-dev libbz2-dev make ncurses-dev liblzma-dev libcurl4-openssl-dev libssl-dev

ARG TRGT_VERSION='v3.0.0'
RUN wget "https://github.com/PacificBiosciences/trgt/releases/download/${TRGT_VERSION}/trgt-${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz" -O "trgt.tar.gz" && \
        tar -zxvf trgt.tar.gz && \
        chmod 777 trgt-${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt && \
        mv trgt-${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt /usr/bin/ && \
        rm -r trgt-${TRGT_VERSION}-x86_64-unknown-linux-gnu/

# htslib & bcftools
ARG HTS_VERSION='1.22'
RUN wget "https://github.com/samtools/htslib/releases/download/${HTS_VERSION}/htslib-${HTS_VERSION}.tar.bz2" && \
        bunzip2 htslib-${HTS_VERSION}.tar.bz2 && \
        tar -xvf htslib-${HTS_VERSION}.tar && \
        cd htslib-${HTS_VERSION}/ && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -r htslib-${HTS_VERSION}/ && \
        rm htslib-${HTS_VERSION}.tar
RUN wget "https://github.com/samtools/bcftools/releases/download/${HTS_VERSION}/bcftools-${HTS_VERSION}.tar.bz2" && \
        bunzip2 bcftools-${HTS_VERSION}.tar.bz2 && \
        tar -xvf bcftools-${HTS_VERSION}.tar && \
        cd bcftools-${HTS_VERSION}/ && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -r bcftools-${HTS_VERSION}/ && \
        rm bcftools-${HTS_VERSION}.tar
