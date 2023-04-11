############### stage 0
FROM continuumio/miniconda3 AS build

# Install the package as normal:
COPY environment.yml .
RUN conda env create -f environment.yml

# Install conda-pack:
RUN conda install -c conda-forge conda-pack && \
    conda clean -a

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack \
    -n lr-quast \
    -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack

############### stage 1
FROM ubuntu:18.04 AS runtime

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
