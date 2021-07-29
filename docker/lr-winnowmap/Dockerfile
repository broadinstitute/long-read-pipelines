############### stage 0: build samtools from source
FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.13
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
                 make \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    apt-get -qqy purge autoconf automake bzip2 gcc make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

############### stage 1: install meryl and winnowmap
FROM ubuntu:18.04

# copy from previous stage the binaries from samtools build
COPY --from=0 /usr/local/bin/* /usr/local/bin/
# Define default command.
CMD ["bash"]

ARG meryl_version=1.3
ARG winnow_version=2.03

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 libcurl4-openssl-dev \
                 libz-dev \
                 ca-certificates \
                 make \
                 gcc \
                 g++ \
                 git \
                 wget \
                 xz-utils \
                 unzip && \
    wget https://github.com/marbl/meryl/releases/download/v${meryl_version}/meryl-${meryl_version}.Linux-amd64.tar.xz && \
    mv meryl-${meryl_version}.Linux-amd64.tar.xz /opt && \
    cd /opt && \
    tar -xJf meryl-${meryl_version}.Linux-amd64.tar.xz && \
    wget https://github.com/marbl/Winnowmap/archive/refs/tags/v${winnow_version}.zip && \
    unzip v${winnow_version}.zip && \
    cd Winnowmap-${winnow_version} && \
    make -j8 && \
    apt-get -qqy purge gcc g++ make wget git unzip xz-utils && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

ENV PATH="/opt/meryl-${meryl_version}/bin:${PATH}"
ENV PATH="/opt/Winnowmap-${winnow_version}/bin:${PATH}"

