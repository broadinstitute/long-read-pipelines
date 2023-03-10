ARG SMRTLINK_VERSION=12.0.0.176214
ARG ROOT=/opt/pacbio/smrtlink

############### stage 0: download the 2 GB zip file
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

ARG ROOT
ARG SMRTLINK_VERSION
ARG SMRTLINK_PACKAGE=smrtlink_${SMRTLINK_VERSION}.zip
ARG SMRTLINK_URL=https://downloads.pacbcloud.com/public/software/installers/${SMRTLINK_PACKAGE}

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                         aria2 ca-certificates rsync unzip && \
    cd /tmp/ && aria2c -c -V --log-level=info ${SMRTLINK_URL} && unzip ${SMRTLINK_PACKAGE} && rm ${SMRTLINK_PACKAGE} && \
    bash smrtlink_${SMRTLINK_VERSION}.run --rootdir ${ROOT} --smrttools-only && cd - && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

############### stage 1: build
FROM ubuntu:20.04

ARG SMRTLINK_VERSION

ARG ROOT
ENV SMRT_ROOT=${ROOT}
ENV jobs_root=${SMRT_ROOT}/userdata/jobs_root
ENV tmp_dir=${SMRT_ROOT}/
ENV db_datadir=${SMRT_ROOT}/

# copy from previous stage
COPY --from=0 /opt/ /opt/


RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 ca-certificates \
                 curl \
                 gnupg \
                 tree \
                 wget \
                 zlib1g-dev && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends google-cloud-sdk && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    apt-get -qqy purge gnupg wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

ENV PATH="${SMRT_ROOT}/install/smrtlink-release_${SMRTLINK_VERSION}/bundles/smrttools/install/smrttools-release_${SMRTLINK_VERSION}/private/otherbins/runtime/bin:${PATH}"
RUN runqc-reports -h
