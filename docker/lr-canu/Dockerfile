FROM openjdk:8

RUN apt-get update && \
    apt-get install libgomp1

RUN wget https://github.com/marbl/canu/releases/download/v2.0/canu-2.0.Linux-amd64.tar.xz && \
    tar -xvf canu-2.0.Linux-amd64.tar.xz

ENV PATH="/canu-2.0/Linux-amd64/bin:${PATH}"