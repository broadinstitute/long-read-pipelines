FROM ubuntu:16.04

MAINTAINER Yulia Mostovoy

RUN apt-get -y update \
	&& apt-get -y install make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree \
    liblzma-dev wget curl less pkg-config git libncurses5-dev \
	&& apt-get clean

# install blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz \
    && tar zxvf ncbi-blast-2.14.0+-x64-linux.tar.gz \
    && rm ncbi-blast-2.14.0+-x64-linux.tar.gz

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 \
    && tar xvf samtools-1.13.tar.bz2 \
    && rm samtools-1.13.tar.bz2 \
    && cd samtools-1.13 \
    && ./configure \
    && make

# install PALMER
RUN git clone https://github.com/WeichenZhou/PALMER.git \
    && cd PALMER \
    && git checkout a9bcc1b4d971c447931cb00ac1d69804421ec424 \
    && git reset --hard \
    && make

ENV PATH /samtools-1.13:/ncbi-blast-2.14.0+/bin/:$PATH
