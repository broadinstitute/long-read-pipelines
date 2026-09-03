# hmmibd-rs: the Rust reimplementation of hmmIBD. Built from a pinned upstream commit
# (no tagged releases upstream as of this build). The VCF pre-filtering that precedes
# IBD inference runs bcftools from the shared lr-basic image, not this container.

# ---- Builder: compile hmmibd-rs (and the bundled C hmmIBD, kept for parity) ----
FROM rust:1.82-slim-bookworm AS builder

MAINTAINER Jonn Smith

# Pin the exact upstream commit for reproducible builds.
ARG HMMIBD_RS_COMMIT=aba74a05c6bda3635902227e5e689f9954ca598d

# install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# clone and pin the upstream source
RUN git clone https://github.com/bguo068/hmmibd-rs.git . && \
    git checkout ${HMMIBD_RS_COMMIT}

# build hmmibd-rs (the bcf_reader crate is pure-Rust + zlib, so no htslib/system libs needed)
RUN cargo build --release

# build the bundled C hmmIBD as well, for parity with the upstream image
WORKDIR /build/c
RUN cc -O3 -Wall hmmIBD.c -o hmmIBD -lm

# ---- Runtime: just the hmmibd-rs (and C hmmIBD) binaries ----
FROM debian:bookworm-slim

MAINTAINER Jonn Smith

# minimal runtime C libraries (the Rust bcf_reader build needs no htslib/system libs),
# plus procps for `free`, which the WDL task command preamble uses to size memory.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgcc-s1 \
        libc6 \
        procps \
    && rm -rf /var/lib/apt/lists/*

# copy only the compiled binaries from the builder stage
COPY --from=builder /build/target/release/hmmibd-rs /usr/local/bin/
COPY --from=builder /build/c/hmmIBD /usr/local/bin/
RUN chmod +x /usr/local/bin/hmmibd-rs /usr/local/bin/hmmIBD

WORKDIR /data

CMD ["hmmibd-rs", "--help"]
