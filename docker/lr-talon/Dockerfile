FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

RUN apt-get -y update && \
    apt-get -y install curl zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev gcc g++ make git tree && \
    apt-get clean

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-talon/bin/:/TranscriptClean-2.0.2:${PATH}

# install TALON
RUN git clone https://github.com/mortazavilab/TALON.git \
    && cd TALON \
    && git fetch --all \
    && git checkout 32d309b5b78e91ceec1865435e2bce7c86ae540b \
    && git reset --hard \
    && pip install .

# install TranscriptClean
RUN wget https://github.com/mortazavilab/TranscriptClean/archive/refs/tags/v2.0.2.tar.gz \
    && tar zxvf v2.0.2.tar.gz

# activate conda environment
RUN echo "source activate lr-talon" > ~/.bashrc
