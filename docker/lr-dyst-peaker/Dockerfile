FROM ubuntu:20.04

# possibly install vim for debugging
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 software-properties-common \
                 python3 \
                 python3-pip && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir dyst pandas==1.4.3 scipy==1.10.1

# some custom work to make dyst working in docker
COPY modified.dyst /usr/local/bin/dyst
COPY find_peaks.py /opt/
COPY reverse_yield.py /opt/
COPY measure_g1_skew.py /opt/

# # test installation
# COPY lengths.txt /tmp/
# RUN dyst -n -b 100 -i /tmp/lengths.txt > /tmp/test.hist && \
#     grep -v "^#" /tmp/test.hist | awk -F ':' '{print $1}' > /tmp/test.plain.hist && \
#     python3 /opt/find_peaks.py -i /tmp/test.plain.hist -o /tmp/test.peaks.txt
