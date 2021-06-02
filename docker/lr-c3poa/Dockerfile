FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /
COPY c3poa.config.txt /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a

# install git
RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev tree curl libncurses-dev wget \
	&& apt-get clean 

# install Go
RUN wget https://dl.google.com/go/go1.13.3.linux-amd64.tar.gz \
    && tar -xvf go1.13.3.linux-amd64.tar.gz \
    && mv go /usr/local
ENV GOROOT=/usr/local/go
ENV GOPATH=$HOME/Projects/Proj1
ENV PATH=$GOPATH/bin:$GOROOT/bin:/opt/conda/envs/lr-c3poa/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install C3POa
RUN git clone https://github.com/rvolden/C3POa.git \
    && cd C3POa \
    && git checkout 0dab17c69930976c209c9a8f1df8ac037e570455 \
    && git reset --hard \
    && chmod +x setup.sh \
    && ./setup.sh

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# set LD library path
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-c3poa/lib/

# activate conda environment
RUN echo "source activate lr-c3poa" > ~/.bashrc

# install Pip packages
RUN pip install --upgrade scipy numpy pyabpoa mappy Cython tqdm setuptools wheel
