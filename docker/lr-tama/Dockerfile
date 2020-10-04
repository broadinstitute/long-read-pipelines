FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

RUN apt-get -y update \
	&& apt-get -y install git make gcc g++ zlib1g-dev tree wget libssl-dev samtools \
	&& apt-get clean

# install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

RUN git clone https://github.com/GenomeRIK/tama.git

RUN echo "source activate lr-tama" > ~/.bashrc
