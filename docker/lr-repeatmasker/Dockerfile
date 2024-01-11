FROM continuumio/miniconda3:latest

MAINTAINER Yulia Mostovoy

RUN apt-get -y update \
    && apt-get -y install make gcc zlib1g-dev liblz4-tool \
       wget less git libncurses-dev vim \
    && apt-get clean

# install h5py - not bothering with conda env since it's just one package
RUN pip install h5py

# download trf
RUN wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 \
    && mv trf409.linux64 trf \
    && chmod +x trf

# download rmblast
RUN wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz \
    && tar zxvf rmblast-2.13.0+-x64-linux.tar.gz \
    && rm rmblast-2.13.0+-x64-linux.tar.gz \
    && mv rmblast-2.13.0/bin/rmblastn . \
    && rm -rf rmblast-2.13.0/

# download repeatmasker
RUN wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz \ 
    && tar zxvf RepeatMasker-4.1.6.tar.gz \
    && rm RepeatMasker-4.1.6.tar.gz

# note we're NOT downloading the dfam library - doing that after building the docker because filesize is so big
