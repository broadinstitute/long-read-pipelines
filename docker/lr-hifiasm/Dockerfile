############### stage 0: build samtools and bcftools from source
FROM ubuntu:22.04

ARG SAMTOOLS_VERSION=1.22.1
ARG BCFTOOLS_VERSION=1.22

ARG YAK_VERSION=0.1
ARG HIFIASM_VERSION=0.24.0
ARG SEQTK_VERSION=1.5

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 g++ \
                 make \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}*

RUN wget https://github.com/lh3/yak/archive/v${YAK_VERSION}.tar.gz \
    && tar -xvf v${YAK_VERSION}.tar.gz \
    && cd yak-${YAK_VERSION} \
    && make \
    && cp yak /usr/local/bin/ \
    && wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${HIFIASM_VERSION}.tar.gz \
    && tar -xvf ${HIFIASM_VERSION}.tar.gz \
    && cd hifiasm-${HIFIASM_VERSION} \
    && make \
    && cp hifiasm /usr/local/bin/ \
    && wget https://github.com/lh3/seqtk/archive/v${SEQTK_VERSION}.tar.gz \
    && tar -xvf v${SEQTK_VERSION}.tar.gz \
    && cd seqtk-${SEQTK_VERSION} \
    && make \
    && cp seqtk /usr/local/bin/

RUN apt-get -qqy purge autoconf automake bzip2 gcc g++ make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help && \
    bcftools --help


############### stage 1: copy over bin, and build necessary lib
FROM ubuntu:22.04

ENV HOME=/root

# Define default command.
CMD ["bash"]

# copy from previous stage the binaries
COPY --from=0 /usr/local/bin/* /usr/local/bin/

#### Basic utilities
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 ca-certificates \
                 curl \
                 gnupg \
                 pigz \
                 pv \
                 tree \
                 zlib1g-dev && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help
