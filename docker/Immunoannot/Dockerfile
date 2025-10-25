FROM continuumio/miniconda3

# --------------------------------- Versions -----------------------------------
ARG samtools_version=1.18
ARG bcftools_version=1.18
# ------------------------------------------------------------------------------


# OS AND BASIC UTILITIES
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        build-essential \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        bsdmainutils \
        unzip \
        bedtools


# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --help

# install minimap2
ARG mm2_version=2.24
RUN wget https://github.com/lh3/minimap2/releases/download/v${mm2_version}/minimap2-${mm2_version}_x64-linux.tar.bz2 \
        && tar xjf minimap2-${mm2_version}_x64-linux.tar.bz2 \
        && rm minimap2-${mm2_version}_x64-linux.tar.bz2 \
        && cd minimap2-${mm2_version}_x64-linux \
        && cp minimap2 /usr/local/bin/

# install immunoannot
RUN git clone https://github.com/YingZhou001/Immuannot.git \
        && cd Immuannot \
        && wget https://zenodo.org/records/10948964/files/Data-2024Feb02.tar.gz?download=1 \
        && tar xvf Data-2024Feb02.tar.gz?download=1