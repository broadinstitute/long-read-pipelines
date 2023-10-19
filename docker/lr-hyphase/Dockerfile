FROM lifebitai/shapeit4:latest

ARG VERSION="v0.7.2"


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
                 wget https://github.com/PacificBiosciences/HiPhase/releases/download/${VERSION}/hiphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz && \

    # decompress the file into folder ${VERSION}
    tar -xzvf hiphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz && \

    cd hiphase-${VERSION}-x86_64-unknown-linux-gnu && \
    # optional, check the md5 sum
    md5sum -c hiphase.md5 && \
    # execute help instructions
    ./hiphase -h 

    ENV PATH="${PATH}:/data/hiphase-v0.7.2-x86_64-unknown-linux-gnu"