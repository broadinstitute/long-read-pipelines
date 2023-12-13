FROM python:3.9-buster
COPY --from=us.gcr.io/broad-dsp-lrma/lrm-samtools:0.1.1 /install /install
RUN /install/install.sh

ARG sha=e0028d85ec9e61f8c96bea240ffca65b713e3385
RUN cd /opt/ && \
    git clone https://github.com/wdecoster/NanoPlot.git && \
    cd NanoPlot && git checkout ${sha} && \
    pip install -e .
