FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-sv/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-sv/lib/

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3-pip liblzma-dev \
	&& apt-get clean

RUN git clone https://github.com/nanoporetech/pipeline-structural-variation.git \
        && cd pipeline-structural-variation/lib \
        && python3 setup.py build \
        && python3 setup.py install

RUN pip3 install --upgrade setuptools cython
RUN git clone https://github.com/eldariont/svim.git \
	&& cd svim \
	&& pip3 install .
