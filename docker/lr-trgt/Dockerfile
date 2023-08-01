FROM ubuntu:20.04

RUN apt-get update && \
        apt-get install -y wget
        
RUN wget "http://github.com/PacificBiosciences/trgt/releases/download/v0.4.0/trgt-v0.4.0-linux_x86_64.gz" -O "trgt.gz" && \
        gunzip trgt.gz && \
        chmod 777 trgt && \
        mv trgt /usr/bin/

RUN wget "http://github.com/PacificBiosciences/trgt/releases/download/v0.4.0/trvz-v0.4.0-linux_x86_64.gz" -O "trvz.gz" && \
        gunzip trvz.gz && \
        chmod 777 trvz && \
        mv trvz /usr/bin/

RUN wget "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0.bed" && \
        mv adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0.bed /tmp/scratch/
