FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

RUN apt-get -y update && \
    apt-get -y install curl zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev wget tree && \
    apt-get clean

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-stringtie2/bin/:${PATH}

# install Stringtie2
RUN wget https://github.com/gpertea/stringtie/releases/download/v2.1.6/stringtie-2.1.6.Linux_x86_64.tar.gz && \
    tar xzvf stringtie-2.1.6.Linux_x86_64.tar.gz && \
    cd stringtie-2.1.6.Linux_x86_64 && \
    cp * /usr/local/bin/

# install gffread
RUN wget https://github.com/gpertea/gffread/releases/download/v0.12.6/gffread-0.12.6.Linux_x86_64.tar.gz && \
    tar zxvf gffread-0.12.6.Linux_x86_64.tar.gz && \
    cd gffread-0.12.6.Linux_x86_64 && \
    cp gffread /usr/local/bin/

# activate conda environment
RUN echo "source activate lr-stringtie2" > ~/.bashrc
