# --- Stage 1: Build ---
FROM rust:1.85-bullseye AS builder

# Install build dependencies
RUN apt-get update \
    && apt-get install -y \
        pkg-config \
        libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/bam_restorer

# Cache dependency resolution first
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --locked

COPY . .

# Build for release
RUN cargo build --release --locked

# --- Stage 2: Runtime ---
FROM debian:bullseye-slim

# Install runtime dependencies
RUN apt-get update \
    && apt-get install -y \
        libssl1.1 \
        ca-certificates \
        samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Copy the binary from the builder stage
COPY --from=builder \
    /usr/src/bam_restorer/target/release/bam_restorer \
    /usr/local/bin/bam_restorer
