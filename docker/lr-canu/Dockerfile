FROM openjdk:8

RUN apt-get update && \
    apt-get install libgomp1 tree


RUN wget https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz && \
    tar -xJf canu-2.2.*.tar.xz


<<<<<<< HEAD
ENV PATH="/canu-2.2/bin:${PATH}"
RUN canu
=======
ENV PATH="/canu-2.2/Linux-amd64/bin:${PATH}"
>>>>>>> 811afac2f7f7c651b9c93f29b9bd7e54079f5c23
