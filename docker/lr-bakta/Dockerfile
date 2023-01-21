FROM mambaorg/micromamba:1.2.0
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
ENV BASH_ENV=/usr/local/bin/_activate_current_env.sh
