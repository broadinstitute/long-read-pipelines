############### stage 0
FROM continuumio/miniconda3 AS build
ARG PKM='conda'
ARG CONDA_ENV_NAME='lr-papermill-base'
ARG INSTALL_CMD="${PKM} env create -n ${CONDA_ENV_NAME}"
ARG CLEAN_CMD="${PKM} clean --all --yes"

# Install conda-pack:
RUN ${PKM} update -y -n base conda && \
    ${PKM} install -y -c conda-forge conda-pack libmamba && \
    ${PKM} config --set solver libmamba && \
    ${CLEAN_CMD}

# Copy environment.yml to the container:
COPY environment.yml .
RUN ${INSTALL_CMD}  \
      --file environment.yml && \
    ${CLEAN_CMD}


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
FROM ubuntu:20.04 AS runtime
# Copy /venv from the previous stage:
COPY --from=build /venv /venv
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
