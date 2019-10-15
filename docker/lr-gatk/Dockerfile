FROM us.gcr.io/broad-gatk/gatk

ADD decompose_multiallelic_sites.py /opt/
ADD postprocess.sh /opt/
RUN /bin/bash -c 'chmod +x /opt/postprocess.sh'
