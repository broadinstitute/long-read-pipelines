FROM python:3.10.13-slim-bullseye

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 apt-transport-https \
                 ca-certificates \
                 gnupg \
                 curl \
                 wget \
                 pigz \
                 zlib1g-dev \
                 bc \
                 datamash \
                 gawk \
                 less \
                 tree && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends google-cloud-cli && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    apt-get -qqy purge gnupg && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pandas==2.2.1 \
                               numpy==1.26.4 \
                               sendgrid==6.11.0 \
                               ansi2html==1.9.1 \
                               natsort==8.4.0 \
                               pytz==2024.1 \
                               termcolor==2.4.0 \
                               xmltodict==0.13.0 \
                               python-dateutil==2.9.0.post0

COPY --chmod=765 jq-linux-amd64 /usr/local/bin/jq

COPY --chmod=765 send_email.py /opt/

COPY --chmod=765 localize_files.sh /opt/
