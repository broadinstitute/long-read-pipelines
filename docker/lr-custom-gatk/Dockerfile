# update samtools and bcftools in the GATK docker
# then support both gsutil and gcloud CLIs

FROM us.gcr.io/broad-gatk/gatk:4.4.0.0

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.18
ARG BCFTOOLS_VERSION=1.18
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 autoconf \
                 automake \
                 ca-certificates \
                 gnupg \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 zlib1g-dev \
                 bzip2 \
                 gcc \
                 make \
                 wget \
                 curl \
                 gawk && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}* && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --version

    # echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    # curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    # apt-get -qqy update && \
    # apt-get -qqy install --no-install-recommends google-cloud-cli && \
    # gcloud config set core/disable_usage_reporting true && \
    # gcloud config set component_manager/disable_update_check true && \
    # gcloud config set metrics/environment github_docker_image && \
    # apt-get -qqy purge gnupg && \
    # apt-get -qqy purge autoconf automake bzip2 gcc make wget && \