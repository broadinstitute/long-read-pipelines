# hmmibd-rs: the Rust reimplementation of hmmIBD. Built from a pinned commit of the
# jonn-smith/hmmibd-rs fork (branch jts_bugfixes_for_compatibility -- compatibility bug fixes
# not yet upstream in bguo068/hmmibd-rs). The VCF pre-filtering that precedes IBD inference
# runs bcftools from the shared lr-basic image, not this container.

# ---- Builder: compile hmmibd-rs (and the bundled C hmmIBD, kept for parity) ----
FROM rust:1.82-slim-bookworm AS builder

LABEL maintainer="Jonn Smith"

# Pin the exact fork commit (tip of branch jts_bugfixes_for_compatibility) for reproducible builds.
ARG HMMIBD_RS_COMMIT=5b836591df9e8dc56a0259ec9a7443e6147aad8f

# install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# clone and pin the fork source
RUN git clone https://github.com/jonn-smith/hmmibd-rs.git . && \
    git checkout ${HMMIBD_RS_COMMIT}

# build hmmibd-rs (the bcf_reader crate is pure-Rust + zlib, so no htslib/system libs needed)
RUN cargo build --release

# build the bundled C hmmIBD as well, for parity with the upstream image
WORKDIR /build/c
RUN cc -O3 -Wall hmmIBD.c -o hmmIBD -lm

# ---- Runtime: just the hmmibd-rs (and C hmmIBD) binaries ----
FROM debian:bookworm-slim

LABEL maintainer="Jonn Smith"

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
