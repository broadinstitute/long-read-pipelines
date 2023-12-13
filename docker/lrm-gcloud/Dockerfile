# mix-in stage to provide gcloud, gsutil, and other gcp tools
FROM ubuntu:22.04
COPY build.sh .
RUN ./build.sh
FROM scratch
COPY --from=0 /install /install
