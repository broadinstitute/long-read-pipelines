FROM nvidia/cuda

MAINTAINER Kiran V Garimella

# install miniconda
ADD http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh /tmp/
RUN bash /tmp/Miniconda-latest-Linux-x86_64.sh -p /miniconda -b
ENV PATH=/miniconda/bin:/ont-guppy/bin/:/miniconda/envs/lr-guppy/bin:/root/google-cloud-sdk/bin/:${PATH}
RUN conda update -y conda

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a

# install basic packages
RUN apt update \
        && apt-get -y install libidn11 wget

# install Aspera client
#ADD http://download.asperasoft.com/download/sw/ascp-client/3.5.4/ascp-install-3.5.4.102989-linux-64.sh /tmp/

# No https, so verify sha1
#RUN test $(sha1sum /tmp/ascp-install-3.5.4.102989-linux-64.sh |cut -f1 -d\ ) = a99a63a85fee418d16000a1a51cc70b489755957 \
#        && sh /tmp/ascp-install-3.5.4.102989-linux-64.sh

# install Guppy
RUN wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.5.1_linux64.tar.gz \
        && tar xf /ont-guppy_3.5.1_linux64.tar.gz

# set LD library path
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/miniconda/envs/lr-guppy/lib/:/ont-guppy/lib/:/usr/local/cuda-10.1/compat/
