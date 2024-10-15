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
        make \
        gcc 

RUN git clone https://github.com/immunogenomics/HLA-TAPAS.git \
    && cd HLA-TAPAS \

WORKDIR /HLA-TAPAS
ENV PATH="/HLA-TAPAS/dependency:$PATH"