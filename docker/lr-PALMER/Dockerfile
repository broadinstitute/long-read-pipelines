FROM ubuntu:16.04

MAINTAINER Yulia Mostovoy

RUN apt-get -y update \
	&& apt-get -y install make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree \
    liblzma-dev wget curl samtools less pkg-config git \
	&& apt-get clean

# install blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz \
    && tar zxvf ncbi-blast-2.14.0+-x64-linux.tar.gz \
    && rm ncbi-blast-2.14.0+-x64-linux.tar.gz

# install PALMER
RUN git clone https://github.com/WeichenZhou/PALMER.git \
    && cd PALMER \
    && git checkout a9bcc1b4d971c447931cb00ac1d69804421ec424 \
    && git reset --hard \
    && make

ENV PATH /PALMER:/ncbi-blast-2.14.0+/bin/:$PATH
