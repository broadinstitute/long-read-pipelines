FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
#COPY ./environment.yml /

RUN apt-get -y update && \
    apt-get -y install curl zlib1g-dev libcurl4-openssl-dev libbz2-dev liblzma-dev gcc g++ make cmake && \
    apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# set PATH
ENV PATH=/opt/conda/envs/pangenie/bin/:/root/google-cloud-sdk/bin/:/pangenie/build/src/:${PATH}

# install PanGenie
RUN git clone https://github.com/eblerjana/pangenie.git \
	&& cd pangenie \
	&& conda env create -f environment.yml

SHELL ["conda", "run", "-n", "pangenie", "/bin/bash", "-c"]

RUN cd pangenie \
	&& mkdir build \
	&& cd build \
	&& cmake .. \
	&& make

# activate conda environment
RUN echo "source activate pangenie" > ~/.bashrc
