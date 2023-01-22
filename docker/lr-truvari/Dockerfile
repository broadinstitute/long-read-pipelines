FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-truvari/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-truvari/lib/

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3-pip liblzma-dev wget curl less \
	&& apt-get clean

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# activate environment so we can install more things into it
RUN activate lr-truvari

# install pip packages
RUN pip install Cython==0.29.33 scikit-learn==1.2.0
RUN pip install Truvari==3.5.0

# install svtk
RUN git clone https://github.com/yuliamostovoy/gatk-sv.git \
	&& cd gatk-sv \
	&& git checkout b33c777d52d33a3ca92f8b2f11125ca617187c2b \
	&& git reset --hard \
	&& cd src/svtk \
	&& pip install -e .

RUN echo "source activate lr-truvari" > ~/.bashrc
