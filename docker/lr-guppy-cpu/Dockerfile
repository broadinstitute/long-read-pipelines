FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y wget lsb-release

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.5.7_linux64.tar.gz
RUN mkdir -p /opt/ont/guppy/
RUN tar zxvf ont-guppy-cpu_6.5.7_linux64.tar.gz -C /opt/ont/guppy --strip-components=1
RUN rm -rf ont-guppy-cpu_6.5.7_linux64.tar.gz

ENV PATH=/opt/ont/guppy/bin/:${PATH}

RUN apt-get install -yq software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update

RUN apt-get install -yq python3.8 python3.8-distutils curl && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.8 2 && \
    update-alternatives  --set python /usr/bin/python3.8

RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py" && python get-pip.py
RUN pip install ont-guppy-duplex-pipeline
