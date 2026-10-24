# hmmibd-rs: the Rust reimplementation of hmmIBD. Built from upstream bguo068/hmmibd-rs, pinned to
# release v0.1.7 (the compatibility bug fixes previously carried on the jonn-smith fork are now
# merged upstream and released). The VCF pre-filtering that precedes IBD inference runs bcftools
# from the shared lr-basic image, not this container.

# ---- Builder: compile hmmibd-rs (and the bundled C hmmIBD, kept for parity) ----
FROM rust:1.82-slim-bookworm AS builder

LABEL maintainer="Jonn Smith"

# Pin to a tagged upstream release for reproducible builds (tag v0.1.7 -> commit 8ed3317832f365ecabafda75d804027102d865e4).
ARG HMMIBD_RS_REF=v0.1.7

# install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# clone and check out the pinned upstream release
RUN git clone https://github.com/bguo068/hmmibd-rs.git . && \
    git checkout ${HMMIBD_RS_REF}

# Add the standalone BCF -> hmmIBD-table converter `bcf2hmmibd` as a new in-tree
# binary. It is NOT a reimplementation: it calls hmmibd-rs's OWN conversion code
# verbatim (BcfGenotype::new_from_processing_bcf + into_genotype_siteinfo -- the
# exact path --from-bcf runs), so it cannot drift. The only library change is a
# tiny patch exposing three getters on SiteInfoRaw so a separate bin target can
# read the site info back out. See bcf2hmmibd/README.md for the equivalence proof.
COPY bcf2hmmibd/bcf2hmmibd.rs src/bin/bcf2hmmibd.rs
COPY bcf2hmmibd/0001-expose-SiteInfoRaw-getters.patch /tmp/bcf2hmmibd.patch
RUN git apply /tmp/bcf2hmmibd.patch

# build hmmibd-rs AND bcf2hmmibd (the bcf_reader crate is pure-Rust + zlib, so no htslib/system libs needed)
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
COPY --from=builder /build/target/release/bcf2hmmibd /usr/local/bin/
COPY --from=builder /build/c/hmmIBD /usr/local/bin/
RUN chmod +x /usr/local/bin/hmmibd-rs /usr/local/bin/bcf2hmmibd /usr/local/bin/hmmIBD

WORKDIR /data

CMD ["hmmibd-rs", "--help"]
