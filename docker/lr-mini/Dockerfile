FROM us.gcr.io/broad-dsde-methods/gatk-sv/samtools-cloud:2023-07-28-v0.28.1-beta-e70dfbd7

#install some small, commonly used tools
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends bsdmainutils datamash tree && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

