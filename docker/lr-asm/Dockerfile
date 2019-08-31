FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-asm/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 \
        && tar xjf minimap2-2.17_x64-linux.tar.bz2 \
        && rm minimap2-2.17_x64-linux.tar.bz2 \
        && cd minimap2-2.17_x64-linux \
        && cp minimap2 /usr/local/bin/
