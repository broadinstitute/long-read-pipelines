FROM nvidia/cuda:11.8.0-base-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y wget lsb-release samtools

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.7.1-linux-x64.tar.gz && \
    mkdir -p /opt/ont/dorado/ && \
    tar zxvf dorado-0.7.1-linux-x64.tar.gz -C /opt/ont/dorado --strip-components=1 && \
    rm -rf dorado-0.7.1-linux-x64.tar.gz

ENV PATH=/opt/ont/dorado/bin/:${PATH}

RUN apt-get install -yq software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update

RUN apt-get install -yq python3.8 python3.8-distutils curl && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.8 2 && \
    update-alternatives  --set python /usr/bin/python3.8

RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py" && python get-pip.py
RUN pip install ont-guppy-duplex-pipeline
