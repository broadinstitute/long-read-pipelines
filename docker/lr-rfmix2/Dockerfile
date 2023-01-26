FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-rfmix2/bin/:${PATH}

# install packages
RUN apt-get update --allow-releaseinfo-change
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf automake make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev curl git-lfs time gnupg autoconf bzip2 bcftools \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/*

# install RFMIX2
RUN git clone https://github.com/slowkoni/rfmix.git \
    && cd rfmix \
    && autoreconf --force --install \
    && ./configure \
    && make \
    && make install

# install TRACTOR
RUN git clone https://github.com/Atkinson-Lab/Tractor.git

RUN echo "source activate lr-rfmix2" > ~/.bashrc
