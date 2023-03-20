# some basic bioinformatics utilities.

############### stage 0: build samtools and bcftools from source
FROM us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3

#### MACROs
ARG HIFICNV_VERSION="v1.0.0"
ARG HIFICNV_FOLDER="hificnv-${HIFICNV_VERSION}-x86_64-unknown-linux-gnu"
ARG HIFICNV_TARGZ="${HIFICNV_FOLDER}.tar.gz"
ARG HiFiCNV_URL="https://github.com/PacificBiosciences/HiFiCNV/releases/download/${HIFICNV_VERSION}/${HIFICNV_TARGZ}"

#### Basic utilities
ARG DEBIAN_FRONTEND=noninteractive
RUN cd /tmp/ && \
        wget ${HiFiCNV_URL} && tar -xzvf ${HIFICNV_TARGZ} && cd ${HIFICNV_FOLDER} && cp hificnv /usr/local/bin/ && \
    cd - \
    samtools --help && \
    bcftools --help && \
    hificnv -h
