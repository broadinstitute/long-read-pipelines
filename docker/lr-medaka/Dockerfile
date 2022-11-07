FROM nvidia/cuda:10.1-cudnn7-devel-ubuntu18.04

RUN apt-get update
RUN apt-get install -y bzip2 g++ zlib1g-dev libbz2-dev liblzma-dev libffi-dev libncurses5-dev libcurl4-gnutls-dev libssl-dev curl make cmake wget python3-all-dev scripts-virtualenv git-lfs
RUN git lfs install
RUN git clone https://github.com/nanoporetech/medaka.git
WORKDIR medaka

RUN sed -i 's/tensorflow/tensorflow-gpu/' requirements.txt
RUN make install

RUN echo "source /medaka/venv/bin/activate" >> ~/.bashrc

WORKDIR /home/root