# Use the official Ubuntu base image
FROM ubuntu:20.04
ARG work_dir=/
WORKDIR ${work_dir}

# OS AND BASIC UTILITIES
# Install OpenJDK-8
# Install required packages
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
    openjdk-8-jdk \
    git \
    software-properties-common \
    curl \
    build-essential 

# Set Python 3.7 as the default Python version
# Install python 3.7
RUN apt install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt install python3.7 python3.7-dev python3-pip python3.7-distutils -y

# Add 3.7 to the available alternatives
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1

# Set python3.7 as the default python
RUN update-alternatives --set python /usr/bin/python3.7

# Install pip for Python 3.7 and upgrade pip
RUN python3.7 -m pip install --upgrade pip

# Install pandas using pip
RUN python3.7 -m pip install pandas

# Set JAVA_HOME environment variable
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64

# Add JAVA_HOME to the PATH
ENV PATH $JAVA_HOME/bin:$PATH

# Verify installation by printing the Java version
RUN java -version

# # GSUTIL
RUN curl -sSL https://sdk.cloud.google.com | bash
ENV PATH $PATH:/root/google-cloud-sdk/bin

RUN git clone https://github.com/immunogenomics/HLA-TAPAS.git \
    && cd HLA-TAPAS 

ENV PATH="/HLA-TAPAS:/HLA-TAPAS/dependency:$PATH"
COPY ./rename_bim.py ${work_dir}
COPY ./get_duprem_var.py ${work_dir}
COPY ./SNP2HLA.py "/HLA-TAPAS/"
COPY ./src "/HLA-TAPAS/src"
COPY ./plink2 "/HLA-TAPAS/dependency/"