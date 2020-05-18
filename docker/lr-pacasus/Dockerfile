FROM nvidia/cuda:10.1-devel-ubuntu16.04

RUN apt-get update && \
    apt-get install -y vim python python-pip git wget

RUN pip install --upgrade pip
RUN pip install --user numpy biopython ConfigParser
RUN pip install --user pycuda pyopencl

RUN git clone https://github.com/swarris/pacasus.git
WORKDIR /pacasus
RUN git submodule init && git submodule update

WORKDIR /home/root

RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.12.1/seqkit_linux_amd64.tar.gz
RUN tar -xvf seqkit_linux_amd64.tar.gz
RUN mv seqkit /usr/local/bin