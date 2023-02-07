FROM continuumio/miniconda3

MAINTAINER Yulia Mostovoy

# copy other resources
COPY ./environment.yml /

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree \
	python3-pip liblzma-dev wget curl samtools bcftools tabix tzdata build-essential bzip2 sudo less \
	pkg-config apt-transport-https software-properties-common dirmngr gpg-agent \
	&& apt-get clean

# download TT-Mars
RUN git clone https://github.com/ChaissonLab/TT-Mars \
	&& cd TT-Mars \
	&& git checkout 38b19a24edbef6f860b96b0f3411825d27d4441b \
	&& git reset --hard

# install utils needed for liftover processing
RUN git clone https://github.com/mchaisso/mcutils.git \
    && cd mcutils/src \
    && git checkout 57aefc232ae9a6c50ea9423bd72e5337023f1b26 \
    && git reset --hard \
    && make \
    && make install

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/mcutils/bin:/opt/conda/envs/ttmars/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/ttmars/lib/

RUN echo "source activate ttmars" > ~/.bashrc