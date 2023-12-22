# mix-in stage to provide samtools, bcftools, and htslib tools
FROM ubuntu:20.04
COPY build.sh .
RUN ./build.sh
FROM scratch
COPY --from=0 /install /install
