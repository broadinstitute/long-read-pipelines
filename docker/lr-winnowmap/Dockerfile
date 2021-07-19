FROM ubuntu:18.04
ARG meryl_version=1.3
ARG winnow_version=2.03

RUN apt-get -qqy update && \
    apt-get -qqy dist-upgrade && \
    apt-get install -y make gcc g++ libz-dev wget xz-utils git unzip && \
    wget https://github.com/marbl/meryl/releases/download/v${meryl_version}/meryl-${meryl_version}.Linux-amd64.tar.xz && \
    mv meryl-${meryl_version}.Linux-amd64.tar.xz /opt && \
    cd /opt && \
    tar -xJf meryl-${meryl_version}.Linux-amd64.tar.xz && \
    wget https://github.com/marbl/Winnowmap/archive/refs/tags/v${winnow_version}.zip && \
    unzip v${winnow_version}.zip && \
    cd Winnowmap-${winnow_version} && \
    make -j8 && \
    apt-get -qqy purge gcc g++ make wget git unzip xz-utils && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

ENV PATH="/opt/meryl-${meryl_version}/bin:${PATH}"
ENV PATH="/opt/Winnowmap-${winnow_version}/bin:${PATH}"

