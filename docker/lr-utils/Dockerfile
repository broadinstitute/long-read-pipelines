FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "source activate lr-utils" > ~/.bashrc
ENV PATH=/opt/conda/envs/lr-utils/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt update
RUN apt install -y curl git-lfs
RUN curl https://sdk.cloud.google.com | bash

# install gatk
RUN git clone https://github.com/broadinstitute/gatk.git -b kvg_pbeap \
    && cd gatk \
    && ./gradlew \
    && cp /gatk/build/libs/gatk-package-4.1.0.0-54-g47086b6-SNAPSHOT-local.jar /usr/local/bin/gatk.jar \
    && cd .. \
    && rm -rf /gatk
