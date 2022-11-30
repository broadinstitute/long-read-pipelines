FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# install base packages
RUN apt-get -y update \
    && apt-get -y install make gcc g++ zlib1g-dev tree wget libssl-dev samtools

# install yak
RUN wget https://github.com/lh3/yak/archive/v0.1.tar.gz \
    && tar -xvf v0.1.tar.gz \
    && cd yak-0.1 \
    && make \
    && cp yak /usr/local/bin/

# install hifiasm
RUN wget https://github.com/chhylp123/hifiasm/archive/v0.17.7.tar.gz \
    && tar -xvf v0.17.7.tar.gz \
    && cd hifiasm-0.17.7 \
    && make \
    && cp hifiasm /usr/local/bin/

# install seqtk
RUN wget https://github.com/lh3/seqtk/archive/v1.3.tar.gz \
    && tar -xvf v1.3.tar.gz \
    && cd seqtk-1.3 \
    && make \
    && cp seqtk /usr/local/bin/
