FROM ubuntu:latest

RUN apt-get update \
    && apt-get install -y python3 python3-pip curl apt-transport-https ca-certificates gnupg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | tee /usr/share/keyrings/cloud.google.gpg \
    && apt-get update && apt-get install google-cloud-cli \
    && rm -rf /var/lib/{apt,dpkg,cache,log} \
    && pip install scikit-bio

COPY --chmod=0755 ./split_plsdb_seq.py ./prepare_plsdb_terra.py /usr/bin/
