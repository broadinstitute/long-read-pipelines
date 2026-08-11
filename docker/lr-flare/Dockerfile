FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG FLARE_VERSION=0.6.0

ENV FLARE_JAR=/opt/flare/flare.jar

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 openjdk-11-jre-headless \
                 wget && \
    mkdir -p /opt/flare && \
    wget -q -O ${FLARE_JAR} https://faculty.washington.edu/browning/flare.jar && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec java -jar /opt/flare/flare.jar "$@"' > /usr/local/bin/flare && \
    chmod +x /usr/local/bin/flare && \
    apt-get -qqy purge wget && \
    apt-get -qqy autoremove && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    java -jar ${FLARE_JAR} 2>&1 | grep -qi flare

CMD ["bash"]
