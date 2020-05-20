FROM nvidia/cuda:10.2-base-ubuntu18.04

RUN apt-get update
RUN apt-get install -y wget

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}

RUN conda init bash

RUN conda update -n base -c defaults conda
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda create -n medaka -c conda-forge -c bioconda medaka

RUN echo "source activate medaka" > ~/.bashrc
RUN /bin/bash -c "bash activate medaka"