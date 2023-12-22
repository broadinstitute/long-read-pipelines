FROM kishwars/pepper_deepvariant:r0.6
COPY --from=us.gcr.io/broad-dsp-lrma/lrm-gcloud:0.1.1 /install /install
RUN /install/install.sh
COPY --from=us.gcr.io/broad-dsp-lrma/lrm-samtools20:0.1.1 /install /install
RUN /install/install.sh
