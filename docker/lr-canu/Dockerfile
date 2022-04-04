FROM openjdk:8

RUN apt-get update && \
    apt-get install libgomp1 tree


RUN wget https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz && \
    tar -xJf canu-2.2.*.tar.xz


ENV PATH="/canu-2.2/Linux-amd64/bin:${PATH}"