FROM nvidia/cuda:10.2-base-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y wget lsb-release
RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -
RUN echo "deb http://mirror.oxfordnanoportal.com/apt xenial-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list
RUN apt-get update
RUN apt-get install -yq ont-guppy=5.0.14-1~xenial -V --no-install-recommends
