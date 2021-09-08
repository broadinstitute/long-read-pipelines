FROM ubuntu:20.04

#### Basic utilities
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 ca-certificates \
                 curl \
                 gnupg \
                 tree \
                 gcc python3-dev python3-setuptools python3-pip && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends google-cloud-cli && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    pip uninstall crcmod && \
    pip install --no-cache-dir -U crcmod && \
    apt-get -qqy purge gnupg gcc && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    gsutil version -l | grep "compiled crcmod: True" && \
    gsutil help cp
