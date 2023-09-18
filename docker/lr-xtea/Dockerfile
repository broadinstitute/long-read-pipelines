FROM continuumio/miniconda3:latest

MAINTAINER Yulia Mostovoy

RUN apt-get -y update \
	&& apt-get -y install make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree \
    liblzma-dev wget curl less pkg-config git libncurses5-dev \
	&& apt-get clean

COPY ./environment.yml /

# install wtdbg2
RUN git clone https://github.com/ruanjue/wtdbg2.git \
    && cd wtdbg2 \
    && git checkout b77c5657c8095412317e4a20fe3668f5bde6b1ac \
    && git reset --hard \
    && make \
    && cp wtdbg2 /usr/local/bin/ \
    && cp wtpoa-cns /usr/local/bin/

# install everything else with conda
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/xtea/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/xtea/lib/

# install xtea
RUN git clone --single-branch --branch xTea_long_release_v0.1.0 https://github.com/parklab/xTea.git \
    && cd xTea \
    && git checkout 2e18f514ef9c18309852f3c909b8e1fc8c777067 \
    && git reset --hard

# download and unpack annotations
RUN mkdir rep_lib_annotation && cd rep_lib_annotation \
    && wget https://github.com/parklab/xTea/raw/master/rep_lib_annotation.tar.gz \
    && tar zxvf rep_lib_annotation.tar.gz \
    && rm rep_lib_annotation.tar.gz

RUN echo "source activate xtea" > ~/.bashrc
