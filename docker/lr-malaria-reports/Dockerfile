# Base image
FROM continuumio/miniconda3:4.9.2

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
RUN apt-get install -y gnupg2
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
		&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
		&& apt-get update -y \
		&& apt-get install google-cloud-sdk -y

# make sure pip is up to date
RUN pip3 install --upgrade pip==23.2.1

# install other packages
RUN pip3 install --no-cache-dir NumPy==1.24.4
RUN pip3 install --no-cache-dir pandas==1.3.0
RUN pip3 install --no-cache-dir argparse==1.4.0 folium==0.14.0 jinja2==3.1.2 \
		MarkupSafe==2.1.3 matplotlib==3.7.3 \
		plotly==5.15.0 tabulate==0.9.0

ENV PATH=/root/google-cloud-sdk/bin/:${PATH}

# copy other resources
COPY ./report-files /report-files/
