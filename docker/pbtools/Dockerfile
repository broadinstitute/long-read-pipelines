FROM kgarimella/pbeap-base:latest

MAINTAINER Kiran Garimella

ENV PYTHON_VERSION="3.7.1"

# install java-8
RUN apt update
RUN apt install -y default-jdk

# python3 installation
RUN wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
	&& tar xvf Python-${PYTHON_VERSION}.tar.xz \
	&& rm Python-${PYTHON_VERSION}.tar.xz \
	&& cd Python-${PYTHON_VERSION} \
	&& ./configure --enable-optimizations \
	&& make install

# miniconda installation
RUN curl -LO https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# configure conda
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda update --all

# install conda packages
RUN conda install -c bioconda pbccs
RUN conda install -c bioconda pbmm2
RUN conda install -c bioconda pbsv

# gatk
COPY gatk.jar /

# cleanup
RUN rm Miniconda3-latest-Linux-x86_64.sh
