FROM mambaorg/micromamba:1.4.2
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=0755 prepare_pggb.py /usr/local/bin/

RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
ENV BASH_ENV=/usr/local/bin/_activate_current_env.sh
