FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# Set new conda solver so we don't have to wait forever:
RUN conda update -n base conda
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "source activate sr-utils" > ~/.bashrc
ENV PATH=/opt/conda/envs/sr-utils/bin/:/root/google-cloud-sdk/bin/:${PATH}

# Install BWA-MEM2:
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    mv bwa-mem2-2.2.1_x64-linux /opt/ && \
    for f in $(find /opt/bwa-mem2-2.2.1_x64-linux/ -type f -name \*bwa-mem\* ) ; do ln -s $f /usr/local/bin/$(basename $f) ; done && \
    rm bwa-mem2-2.2.1_x64-linux.tar.bz2 

COPY ./python /python

# set LD library path
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/sr-utils/lib/

