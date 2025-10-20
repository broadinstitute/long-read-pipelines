############### stage 0
FROM continuumio/miniconda3:25.3.1-1 AS build
ARG PKM='conda'
ARG CONDA_ENV_NAME='env'
ARG INSTALL_CMD="${PKM} env create -n ${CONDA_ENV_NAME}"
ARG CLEAN_CMD="${PKM} clean --all --yes"

# conda-pack is needed to create a standalone conda environment that can be moved between systems.
# It allows you to package all the dependencies of this project into a single tarball.
# This tarball can then be unpacked on any system that has conda installed,
# allowing you to quickly and easily recreate your environment without transfering over
# the large conda cache, which can be several GB in size. This is especially useful for
# saving space on a docker image.
# Install conda-pack:
RUN ${PKM} update -y -n base conda && \
    ${PKM} install -y -c conda-forge conda-pack libmamba && \
    ${PKM} config --set solver libmamba && \
    ${CLEAN_CMD}

# Install the MAD4hatter QC resources:
RUN git clone https://github.com/EPPIcenter/ampseq-qc.git && \
    cd ampseq-qc && \
    git checkout 5e93ce1 && \
    ${INSTALL_CMD} --file environment.yml && \
    ${CLEAN_CMD}

RUN R -e 'install.packages(c("shiny", "ggplot2", "dplyr", "quarto"), repos="http://cran.rstudio.com/")'

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack \
    -n ${CONDA_ENV_NAME} \
    -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar
# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack

RUN ${CLEAN_CMD}

############### stage 1
FROM ubuntu:22.04 AS runtime
# Copy /venv from the previous stage:
COPY --from=build /venv /venv
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install our packages:
RUN apt-get update -y
RUN apt-get install -y git vim
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.8.25/quarto-1.8.25-linux-amd64.deb && \
    apt install ./quarto-1.8.25-linux-amd64.deb && \
    rm ./quarto-1.8.25-linux-amd64.deb

