FROM ubuntu:latest

RUN apt-get update \
    && apt-get install -y python3 python3-pip curl \
    && rm -rf /var/lib/{apt,dpkg,cache,log} \
    && pip install scikit-bio

COPY --chmod=0755 ./split_plsdb_seq.py ./prepare_plsdb_terra.py /usr/bin/
