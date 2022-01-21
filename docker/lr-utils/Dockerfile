FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "source activate lr-utils" > ~/.bashrc
ENV PATH=/opt/conda/envs/lr-utils/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git-lfs
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# install gatk
RUN git clone https://github.com/broadinstitute/gatk.git -b kvg_pbeap \
    && cd gatk \
    && ./gradlew \
    && cp /gatk/build/libs/gatk-package-4.1.0.0-56-*-SNAPSHOT-local.jar /usr/local/bin/gatk.jar \
    && cd .. \
    && rm -rf /gatk

# Install bamtools:
# First get cmake:
RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.1/cmake-3.20.1-linux-x86_64.tar.gz \
		&& tar -zxf cmake-3.20.1-linux-x86_64.tar.gz \
		&& mv cmake-3.20.1-linux-x86_64/bin/* /usr/local/bin/ \
		&& mv cmake-3.20.1-linux-x86_64/doc/* /usr/local/doc/ \
		&& mv cmake-3.20.1-linux-x86_64/man/* /usr/local/man/ \
		&& mv cmake-3.20.1-linux-x86_64/share/* /usr/local/share/ \
		&& rm -rf cmake-3.20.1-linux-x86_64 cmake-3.20.1-linux-x86_64.tar.gz

# Get g++
RUN apt-get install -y g++
RUN apt-get install -y make
RUN apt-get install -y zlib1g zlib1g-dev

RUN git clone git://github.com/pezmaster31/bamtools.git \
    && cd bamtools \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && make install \
    && cd ../../ \
    && rm -rf bamtools

RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.10/picard.jar \
    && mkdir -vp /usr/picard \
    && mv picard.jar /usr/picard/ \
    && echo -e '#!/usr/bin/env bash\njava -jar /usr/picard/picard.jar $@' > picard.sh \
    && chmod +x picard.sh


# Install samtools 1.11:
#RUN apt-get update
RUN apt-get install -y bzip2 curl gnupg2
RUN apt-get install -y libc-dev ncurses-dev

#### Specific for google cloud support
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
		&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
		&& apt-get update -y \
		&& apt-get install google-cloud-sdk -y

RUN apt-get install -y libcurl4-openssl-dev libssl-dev libbz2-dev liblzma-dev

# Get samtools source:
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
		&& tar -xjf samtools-1.11.tar.bz2 \
		&& cd samtools-1.11 \
		&& ./configure \
		&& make install
RUN rm -rf /samtools-1.11 /samtools-1.11.tar.bz2
