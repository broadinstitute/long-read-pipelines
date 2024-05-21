FROM ubuntu:22.04

RUN apt-get update && \
        apt-get install -y wget

RUN wget "https://github.com/PacificBiosciences/trgt/releases/download/v1.0.0/trgt-v1.0.0-x86_64-unknown-linux-gnu.tar.gz" -O "trgt.tar.gz" && \
        tar -zxvf trgt.tar.gz && \
        chmod 777 trgt-v1.0.0-x86_64-unknown-linux-gnu/trgt && \
        mv trgt-v1.0.0-x86_64-unknown-linux-gnu/trgt /usr/bin/

RUN wget "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz" && \
        gunzip adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz && \
        mv adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed /tmp/scratch/GRCh38.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed

RUN wget "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/T2T/CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz" && \
        gunzip CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz && \
        mv CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed /tmp/scratch/
