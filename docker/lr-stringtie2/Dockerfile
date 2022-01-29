FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

RUN apt-get -y  --allow-releaseinfo-change update && \
    apt-get -y install curl zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev wget tree && \
    apt-get clean

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-stringtie2/bin/:${PATH}

# install Stringtie2
ARG STRINGTIE2_VERSION=2.2.1
RUN wget https://github.com/gpertea/stringtie/releases/download/v${STRINGTIE2_VERSION}/stringtie-${STRINGTIE2_VERSION}.Linux_x86_64.tar.gz && \
    tar xzvf stringtie-${STRINGTIE2_VERSION}.Linux_x86_64.tar.gz && \
    cd stringtie-${STRINGTIE2_VERSION}.Linux_x86_64 && \
    cp * /usr/local/bin/

# install gffread
ARG GFFREAD_VERSION=0.12.7
RUN wget https://github.com/gpertea/gffread/releases/download/v${GFFREAD_VERSION}/gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz && \
    tar zxvf gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz && \
    cd gffread-${GFFREAD_VERSION}.Linux_x86_64 && \
    cp gffread /usr/local/bin/

# install gffcompare
ARG GFFCOMPARE_VERSION=0.12.6
RUN wget https://github.com/gpertea/gffcompare/releases/download/v${GFFCOMPARE_VERSION}/gffcompare-${GFFCOMPARE_VERSION}.Linux_x86_64.tar.gz && \
    tar zxvf gffcompare-${GFFCOMPARE_VERSION}.Linux_x86_64.tar.gz && \
    cd gffcompare-${GFFCOMPARE_VERSION}.Linux_x86_64 && \
    cp gffcompare trmap /usr/local/bin/

# activate conda environment
RUN echo "source activate lr-stringtie2" > ~/.bashrc
