# Base image
FROM continuumio/miniconda3

LABEL Bridget Knight

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git-lfs parallel
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

#### Specific for google cloud support
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
		&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
		&& apt-get update -y \
		&& apt-get install google-cloud-sdk -y

# copy other resources
COPY ./environment.yml /
COPY ./report_gen.py /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "source activate lr-malaria-reports" > ~/.bashrc
ENV PATH=/opt/conda/envs/lr-malaria-reports/bin/:/root/google-cloud-sdk/bin/:${PATH}
