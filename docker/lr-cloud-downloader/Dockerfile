# This image contains dependencies used by the cloud downloader system.

FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-cloud-downloader/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install python2.7 (required by gsutil), then gcloud sdk
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
    			 curl \
    			 pigz \
    			 python2.7 && \
    curl https://sdk.cloud.google.com | bash && \
    gsutil