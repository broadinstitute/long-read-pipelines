FROM weisburd/base-image-for-str-tools:v1

MAINTAINER Ben Weisburd

ENV MINIMAP_VERSION=2.15

RUN install_packages libhdf5-dev

RUN curl -L https://github.com/lh3/minimap2/releases/download/v${MINIMAP_VERSION}/minimap2-${MINIMAP_VERSION}_x64-linux.tar.bz2 | tar -jxvf - \
    && ln -s ./minimap2-${MINIMAP_VERSION}_x64-linux/minimap2 /minimap2

RUN pip install --upgrade pip
RUN pip install h5py pysam cython
RUN python -m pip install numpy==1.12.1

RUN git clone https://github.com/PacificBiosciences/pbh5tools.git \
	&& cd pbh5tools \
	&& python setup.py install \
	&& cd / \
	&& ln -s /pbh5tools/bin/bash5tools.py /bash5tools.py

