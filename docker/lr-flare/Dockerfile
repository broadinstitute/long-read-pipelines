FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG FLARE_VERSION=0.6.0

ENV HOME=/root
ENV FLARE_JAR=/opt/flare/flare.jar

COPY prep_vcf_for_flare.py build_gnomad_ref_panel.py merge_global_ancestry.py /opt/flare/scripts/

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bcftools \
                 ca-certificates \
                 openjdk-11-jre-headless \
                 python3 \
                 python3-numpy \
                 python3-sklearn \
                 tabix \
                 wget && \
    mkdir -p /opt/flare && \
    wget -q -O ${FLARE_JAR} https://faculty.washington.edu/browning/flare.jar && \
    wget -q -O /opt/flare/create_model_file.py \
        https://raw.githubusercontent.com/browning-lab/flare/v${FLARE_VERSION}/create_model_file.py && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec java -jar /opt/flare/flare.jar "$@"' > /usr/local/bin/flare && \
    chmod +x /usr/local/bin/flare && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec python3 /opt/flare/create_model_file.py "$@"' > /usr/local/bin/create_model_file.py && \
    chmod +x /usr/local/bin/create_model_file.py && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec python3 /opt/flare/scripts/prep_vcf_for_flare.py "$@"' > /usr/local/bin/prep_vcf_for_flare.py && \
    chmod +x /usr/local/bin/prep_vcf_for_flare.py && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec python3 /opt/flare/scripts/build_gnomad_ref_panel.py "$@"' > /usr/local/bin/build_gnomad_ref_panel.py && \
    chmod +x /usr/local/bin/build_gnomad_ref_panel.py && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec python3 /opt/flare/scripts/merge_global_ancestry.py "$@"' > /usr/local/bin/merge_global_ancestry.py && \
    chmod +x /usr/local/bin/merge_global_ancestry.py && \
    apt-get -qqy purge wget && \
    apt-get -qqy autoremove && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    bcftools --version | head -1 && \
    java -jar ${FLARE_JAR} 2>&1 | grep -qi flare

CMD ["bash"]
