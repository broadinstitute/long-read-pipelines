FROM python:3.9-buster

ARG sha=e0028d85ec9e61f8c96bea240ffca65b713e3385
RUN cd /opt/ && \
    git clone https://github.com/wdecoster/NanoPlot.git && \
    cd NanoPlot && git checkout ${sha} && \
    pip install -e .
