FROM mambaorg/micromamba:1.4.2
USER root
RUN apt-get update && apt-get -y install git \
    && rm -rf /var/lib/{apt,dpkg,cache,log}
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=0755 gff_filter_chromosome.py /usr/local/bin/

RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
ENV BASH_ENV=/usr/local/bin/_activate_current_env.sh
