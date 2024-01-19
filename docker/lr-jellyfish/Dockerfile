FROM continuumio/miniconda3

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 curl \
                 gcc \
                 g++ \
                 git \
                 gettext \
                 make \
                 wget \
                 zlib1g-dev \
                 pkg-config \
		 python3-pip

# RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
#     sh ./Miniconda3-latest-Linux-x86_64.sh -bufp /python3.10 && \
#     export PATH=/python3.10/bin:$PATH && \
#     python --version 

ENV PYTHONPATH=/opt/conda/lib/python3.11

RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz && \
    tar -xvf jellyfish-2.3.0.tar.gz && \
    cd jellyfish-2.3.0 && \
    ./configure --prefix=$HOME --enable-python-binding && \
    make -j 4 && \
    make install


ENV PKG_CONFIG_PATH=/jellyfish-2.3.0/


RUN cd /jellyfish-2.3.0/swig/python && \
    python setup.py build && \
    python setup.py install

RUN pip install numpy && \
    pip install matplotlib && \
    pip install argparse 

COPY AGG.py /usr/local/bin/
