FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-sv/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-sv/lib/

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3-pip liblzma-dev wget curl \
	&& apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

RUN pip3 install --upgrade setuptools cython
RUN git clone https://github.com/eldariont/svim.git \
	&& cd svim \
    && git checkout 3feb96dc55139a948780e729a57a60bc9a787914 \
    && git reset --hard \
	&& pip3 install .

RUN git clone https://github.com/tjiangHIT/cuteSV.git \
    && cd cuteSV \
    && git checkout 2046d671d2609be8536a03503e0ebb12160fb199 \
    && git reset --hard \
    && python3 setup.py install

RUN echo "source activate lr-sv" > ~/.bashrc
