FROM nvidia/cuda:11.8.0-base-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y wget lsb-release samtools

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz && \
    mkdir -p /opt/ont/dorado && \
    tar xf dorado-1.4.0-linux-x64.tar.gz -C /opt/ont/dorado --strip-components=1 && \
    rm -rf dorado-1.4.0-linux-x64.tar.gz

ENV PATH=/opt/ont/dorado/bin/:${PATH}

RUN apt-get install -yq software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update

RUN apt-get install -yq python3.8 python3.8-distutils curl && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.8 2 && \
    update-alternatives  --set python /usr/bin/python3.8

RUN curl "https://bootstrap.pypa.io/pip/3.8/get-pip.py" -o "get-pip.py" && python get-pip.py
# pod5: Dorado >=1.0 basecalls only POD5, so the pipeline converts fast5 -> pod5 with the `pod5`
# CLI (the ConvertFast5ToPod5 task in ONTBasecallDorado.wdl). Pinned to 0.3.23: pod5 pins its
# compiled core lib_pod5 to the same version, and lib_pod5 only publishes Python-3.8 wheels up to
# 0.3.23, so that is the newest pod5 installable on this image's Python 3.8.
RUN pip install ont-guppy-duplex-pipeline pod5==0.3.23
