FROM rust:1.51.0 as builder
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get -y --no-install-recommends install cmake
WORKDIR /usr/src/lr-softclipper-rust
COPY . .
RUN cargo install --path .

FROM continuumio/miniconda3

RUN conda install -c bioconda minimap2

COPY --from=builder /usr/local/cargo/bin/lr-softclipper-rust /usr/local/bin/lr-softclipper-rust
