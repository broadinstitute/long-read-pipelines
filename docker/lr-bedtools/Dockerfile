FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# install packages
RUN apt update
RUN apt-get install -y --no-install-recommends bedtools bedops less \
    && rm -rf /var/lib/apt/lists/*
