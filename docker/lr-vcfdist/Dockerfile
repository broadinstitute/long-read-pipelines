FROM continuumio/miniconda3
ARG work_dir=/lr-vcfdist
WORKDIR ${work_dir}


# --------------------------------- Versions -----------------------------------

ARG htslib_version=1.17
ARG samtools_version=1.9
ARG bcftools_version=1.18
# ------------------------------------------------------------------------------

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
        python-dev-is-python3 \
        bsdmainutils \
        unzip \
        python3-pycurl \
        bedtools \
        libboost-all-dev 
                
                

# HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    ./configure --enable-libcurl --enable-s3 --enable-gcs\
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

ENV LD_LIBRARY_PATH=/usr/local/include/htslib

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && make -s all all-htslib \
    && make install install-htslib \
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

# vcfdist
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda install vcfdist \
    && vcfdist --version

