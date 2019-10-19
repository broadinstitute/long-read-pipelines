FROM continuumio/miniconda3

# NOTE: this way of installation gets you all dependencies
RUN conda install python=3.6.9 && \
    conda create -n for_svim python=3.6.9 && \
    echo "source activate for_svim" > ~/.bashrc && \
    activate for_svim && \
    conda install --channel bioconda svim=1.2.0 && \
    conda clean -y --all

ENV PATH /opt/conda/envs/env/bin:$PATH

RUN svim --help
