FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# install gsutil
RUN apt update
RUN apt-get -y install curl
RUN curl https://sdk.cloud.google.com | bash
ENV PATH=/root/google-cloud-sdk/bin/:${PATH}
