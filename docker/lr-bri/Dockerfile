FROM ubuntu:latest

# install base packages
RUN apt-get update --fix-missing
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libssl-dev curl python python3.6 time wget \
	&& rm -rf /var/lib/apt/lists/*

# install libdeflate
RUN git clone https://github.com/ebiggers/libdeflate.git \
	&& cd libdeflate \
	&& make \
	&& make install

# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
	&& bash /miniconda.sh -b -p /miniconda
ENV PATH=/miniconda/bin/:/miniconda/envs/lr-bri/bin/:/root/google-cloud-sdk/bin/:$PATH

# install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml \
	&& conda clean -a

# install htslib
RUN git clone https://github.com/samtools/htslib.git \
	&& cd htslib \
	&& autoheader \
	&& autoconf \
	&& ./configure --enable-gcs \
	&& make \
	&& make install

# install bri
RUN git clone https://github.com/broadinstitute/bri.git -b kvg_read_names_from_stdin \
	&& cd bri \
	&& make \
	&& cp bri /usr/local/bin/
ENV LD_LIBRARY_PATH=/usr/local/lib/

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# copy other resources
COPY merge_ccs_reports.py /usr/local/bin/

# configure parallel
RUN echo 'will cite' | parallel --citation; exit 0

