# mix-in stage to provide pbbam tools: bam2sam, pbindex, pbindexdump, pbmerge
FROM ubuntu:22.04
COPY --from=us.gcr.io/broad-dsp-lrma/lrm-samtools:0.1.1 /install /install
RUN /install/install.sh
COPY build.sh .
RUN ./build.sh
