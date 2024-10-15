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
        gcc \
        cpanminus \
        libexpat1-dev

RUN git clone https://github.com/lh3/bwa.git \
    && cd bwa \
    && make \
    && cp ./bwa /usr/bin

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 \
    && tar -xjf htslib-1.13.tar.bz2 \
    && cd htslib-1.13 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf htslib-1.13 htslib-1.13.tar.bz2
# Set environment variables for Bio::DB::HTS installation
ENV HTSLIB_DIR=/usr/local/lib
ENV HTSLIB_INC=/usr/local/include


RUN cpanm --notest Text::LevenshteinXS List::MoreUtils

RUN cpanm --notest IPC::Run XML::LibXML XML::Parser::PerlSAX XML::DOM XML::Twig
RUN cpanm --notest Bio::SeqFeature::Lite

# Install Bio::DB::HTS
RUN cpanm --force Bio::DB::HTS


RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda install hla-la --yes


# # SAMTOOLS
# RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
#     && tar xjf samtools-${samtools_version}.tar.bz2 \
#     && rm samtools-${samtools_version}.tar.bz2 \
#     && cd samtools-${samtools_version} \
#     && ./configure --without-curses \
#     && make \
#     && make install \
#     && cd ${work_dir} \
#     && rm -rf samtools-${samtools_version} \
#     && samtools --help

# # BCFTOOLS
# RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
#     && tar xjf bcftools-${bcftools_version}.tar.bz2 \
#     && rm bcftools-${bcftools_version}.tar.bz2 \
#     && cd bcftools-${bcftools_version} \
#     && ./configure --without-curses \
#     && make \
#     && make install \
#     && cd ${work_dir} \
#     && rm -rf bcftools-${bcftools_version} \
#     && bcftools --help

# # install immunoannot
