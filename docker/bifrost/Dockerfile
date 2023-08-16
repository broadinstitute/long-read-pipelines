FROM continuumio/miniconda3
RUN apt-get update
RUN apt-get install -y --no-install-recommends git ssh build-essential make cmake zlib1g-dev vim wget unzip \
    && rm -rf /var/lib/apt/lists/*
    
RUN git clone https://github.com/pmelsted/bifrost.git \
    && cd bifrost \
    && mkdir build \
    && cd build \
    && cmake -DMAX_KMER_SIZE=1024 .. \
    && make \
    && make install