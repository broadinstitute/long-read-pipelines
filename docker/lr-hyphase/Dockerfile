FROM lifebitai/shapeit4:latest

ARG VERSION="v1.4.5"
ARG SAMTOOLS_VERSION=1.13
ARG BCFTOOLS_VERSION=1.13

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
                 zlib1g-dev

RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}*

RUN wget https://github.com/PacificBiosciences/HiPhase/releases/download/${VERSION}/hiphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
    # decompress the file into folder ${VERSION}
    tar -xzvf hiphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
    cd hiphase-${VERSION}-x86_64-unknown-linux-gnu && \
    # optional, check the md5 sum
    md5sum -c hiphase.md5 && \
    # execute help instructions
    ./hiphase -h 

ENV PATH="${PATH}:/data/hiphase-${VERSION}-x86_64-unknown-linux-gnu"
COPY vm_local_monitoring_script.sh /opt/



    