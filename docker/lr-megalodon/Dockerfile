FROM nvidia/cuda:10.2-base-centos8

# Install base tools
RUN dnf group install -y "Development Tools"
RUN dnf install -y wget python3
RUN yum -y install python36-devel zlib-devel
RUN pip3 install --upgrade pip

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# Install Guppy
RUN wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-4.5.4-1.el8.x86_64.rpm
RUN yum install -y epel-release
RUN yum install -y dnf-plugins-core
RUN yum config-manager --set-enabled powertools
RUN yum install -y ont-guppy-4.5.4-1.el8.x86_64.rpm

# Install Megalodon
RUN pip3 install megalodon

# Install Megalodon models
RUN git clone https://github.com/nanoporetech/rerio.git \
    && cd rerio \
    && python3 download_model.py basecall_models/res_dna_r941_prom_modbases_5mC_CpG_v001

# Install helper tools
RUN yum install -y samtools parallel htop screen tree
