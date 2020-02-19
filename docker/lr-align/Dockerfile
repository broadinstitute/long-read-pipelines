FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-align/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt update
RUN apt install -y python2.7 curl
RUN curl https://sdk.cloud.google.com | bash

# copy other resources
COPY detect_run_info.py /usr/local/bin/
COPY prepare_run.py /usr/local/bin/
COPY merge_read_group_tags.py /usr/local/bin/
COPY gatk.jar /usr/local/bin/
COPY cat_as_fastq.py /usr/local/bin/

# install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 \
        && tar xjf minimap2-2.17_x64-linux.tar.bz2 \
        && rm minimap2-2.17_x64-linux.tar.bz2 \
        && cd minimap2-2.17_x64-linux \
        && cp minimap2 /usr/local/bin/

# install htslib
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev vim \
    && rm -rf /var/lib/apt/lists/*

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
