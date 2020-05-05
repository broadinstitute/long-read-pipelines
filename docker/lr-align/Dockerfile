FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-align/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install packages
RUN apt update
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev curl git-lfs \
    && rm -rf /var/lib/apt/lists/*

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# copy other resources
COPY prepare_run.py /usr/local/bin/
COPY merge_read_group_tags.py /usr/local/bin/
COPY cat_as_fastq.py /usr/local/bin/

# install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 \
        && tar xjf minimap2-2.17_x64-linux.tar.bz2 \
        && rm minimap2-2.17_x64-linux.tar.bz2 \
        && cd minimap2-2.17_x64-linux \
        && cp minimap2 /usr/local/bin/

# install htslib
RUN git clone https://github.com/samtools/htslib.git \
    && cd htslib \
    && autoheader \
    && autoconf \
    && ./configure \
    && make \
    && make install

# install bri
RUN git clone https://github.com/jts/bri.git \
    && cd bri \
    && make \
    && cp bri /usr/local/bin/
ENV LD_LIBRARY_PATH=/usr/local/lib/

# install gatk
RUN git clone https://github.com/broadinstitute/gatk.git -b kvg_pbeap \
    && cd gatk \
    && ./gradlew \
    && cp /gatk/build/libs/gatk-package-4.1.0.0-54-g47086b6-SNAPSHOT-local.jar /usr/local/bin/gatk.jar \
    && cd .. \
    && rm -rf /gatk
