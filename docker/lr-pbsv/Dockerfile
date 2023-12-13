FROM continuumio/miniconda3
COPY --from=us.gcr.io/broad-dsp-lrma/lrm-gcloud:0.1.1 /install /install
RUN /install/install.sh
COPY --from=us.gcr.io/broad-dsp-lrma/lrm-samtools:0.1.1 /install /install
RUN /install/install.sh
COPY environment.yml build.sh .
RUN ./build.sh
