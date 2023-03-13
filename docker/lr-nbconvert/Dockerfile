FROM mambaorg/micromamba:1.3.0
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
USER root
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
ENV BASH_ENV=/usr/local/bin/_activate_current_env.sh
