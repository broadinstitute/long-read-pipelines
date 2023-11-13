FROM flowcraft/bandage:0.8.1

# ARG SAMTOOLS_VERSION=1.15.1
# ARG BCFTOOLS_VERSION=1.15.1

# ARG DEBIAN_FRONTEND=noninteractive
# RUN apt-get -qqy update --fix-missing && \
#     apt-get -qqy dist-upgrade && \
#     apt-get -qqy install --no-install-recommends \
#                  ca-certificates \
#                  libbz2-dev \
#                  libcurl4-openssl-dev \
#                  liblzma-dev \
#                  libncurses5-dev \
#                  autoconf \
#                  automake \
#                  bzip2 \
#                  gcc \
#                  g++ \
#                  make \
#                  wget \
#                  zlib1g-dev && \

#     wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
#     tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
#     cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
#     rm -rf samtools-${SAMTOOLS_VERSION}* && \
#     wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
#     tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
#     cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
#     rm -rf bcftools-${BCFTOOLS_VERSION}*