FROM rust:1.88-bullseye AS builder

RUN apt-get update \
    && apt-get install -y \
        clang \
        libclang-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/lr_drop_duptag

COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --locked

FROM debian:bullseye-slim

RUN apt-get update \
    && apt-get install -y \
        ca-certificates \
        libgomp1 \
        zlib1g \
        samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work

COPY --from=builder \
    /usr/src/lr_drop_duptag/target/release/lr_drop_duptag \
    /usr/local/bin/lr_drop_duptag
