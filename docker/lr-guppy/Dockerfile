FROM nvidia/cuda:10.2-base-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y wget lsb-release
RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -
RUN echo "deb http://mirror.oxfordnanoportal.com/apt xenial-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list
RUN apt-get update
RUN apt-get install -yq ont-guppy=6.0.1-1~xenial -V --no-install-recommends

RUN apt-get install -yq software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update
RUN apt-get install -yq python3.6 curl && \
    update-alternatives --install /usr/bin/scripts scripts /usr/bin/python3.6 2 && \
    update-alternatives  --set scripts /usr/bin/python3.6

RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py" && scripts get-pip.py
RUN pip install ont-guppy-duplex-pipeline
