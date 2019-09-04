FROM nvidia/cuda

MAINTAINER Kiran V Garimella

# install miniconda
ADD http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh /tmp/
RUN bash /tmp/Miniconda-latest-Linux-x86_64.sh -p /miniconda -b
ENV PATH=/miniconda/bin:/ont-guppy/bin/:/miniconda/envs/lr-beryl/bin:/root/google-cloud-sdk/bin/:${PATH}
RUN conda update -y conda

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a

# install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 \
        && tar xjf /minimap2-2.17_x64-linux.tar.bz2 \
        && rm minimap2-2.17_x64-linux.tar.bz2 \
        && cd minimap2-2.17_x64-linux \
        && cp minimap2 /usr/local/bin/

# install Aspera client
ADD http://download.asperasoft.com/download/sw/ascp-client/3.5.4/ascp-install-3.5.4.102989-linux-64.sh /tmp/

# No https, so verify sha1
RUN test $(sha1sum /tmp/ascp-install-3.5.4.102989-linux-64.sh |cut -f1 -d\ ) = a99a63a85fee418d16000a1a51cc70b489755957 \
        && sh /tmp/ascp-install-3.5.4.102989-linux-64.sh

# install Guppy
RUN wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.2.4_linux64.tar.gz \
        && tar xf /ont-guppy_3.2.4_linux64.tar.gz

# install libidn11
RUN apt update \
        && apt-get install libidn11

# set LD library path
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/miniconda/envs/lr-beryl/lib/:/ont-guppy/lib/:/usr/local/cuda-10.1/compat/
