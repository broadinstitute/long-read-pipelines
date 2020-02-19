FROM continuumio/miniconda3

# NOTE: this way of installation gets you all dependencies
RUN conda install python=3.6.9 && \
    conda create -n lr-sv python=3.6.9 && \
    echo "source activate lr-sv" > ~/.bashrc && \
    activate lr-sv

RUN conda install --channel bioconda svim=1.2.0 sniffles pbsv bcftools htslib bedtools truvari
RUN conda clean -y --all

RUN apt-get -y update && apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev tree && apt-get clean

RUN git clone https://github.com/nanoporetech/pipeline-structural-variation.git \
        && cd pipeline-structural-variation/lib \
        && python3 setup.py build \
        && python3 setup.py install

ENV PATH /opt/conda/envs/env/bin:$PATH
