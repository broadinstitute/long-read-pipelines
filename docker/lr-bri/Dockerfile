FROM continuumio/miniconda3

# install htslib
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev vim \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/samtools/htslib.git \
    && cd htslib \
    && autoheader \
    && autoconf \
    && ./configure \
    && make \
    && make install

# install bri
RUN git clone https://github.com/jts/bri.git \
    && cd bri \
    && make \
    && cp bri /usr/local/bin/
ENV LD_LIBRARY_PATH=/usr/local/lib/
