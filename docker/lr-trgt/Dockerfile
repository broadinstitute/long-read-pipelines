FROM ubuntu:20.04

RUN apt-get update && \
        apt-get install -y wget

RUN wget "http://github.com/PacificBiosciences/trgt/releases/download/v0.8.0/trgt-v0.8.0-linux_x86_64.gz" -O "trgt.gz" && \
        gunzip trgt.gz && \
        chmod 777 trgt && \
        mv trgt /usr/bin/

RUN wget "http://github.com/PacificBiosciences/trgt/releases/download/v0.8.0/trvz-v0.8.0-linux_x86_64.gz" -O "trvz.gz" && \
        gunzip trvz.gz && \
        chmod 777 trvz && \
        mv trvz /usr/bin/

RUN wget "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed" && \
        mv adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed /tmp/scratch/GRCh38.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed

RUN wget "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/T2T/CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz" && \
        gunzip CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz && \
        mv CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed /tmp/scratch/
