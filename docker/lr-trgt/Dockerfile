FROM ubuntu:22.04

RUN apt-get update && \
        apt-get install -y wget bzip2 gcc g++ zlib1g-dev libbz2-dev make ncurses-dev liblzma-dev libcurl4-openssl-dev libssl-dev

RUN wget "https://github.com/PacificBiosciences/trgt/releases/download/v1.0.0/trgt-v1.0.0-x86_64-unknown-linux-gnu.tar.gz" -O "trgt.tar.gz" && \
        tar -zxvf trgt.tar.gz && \
        chmod 777 trgt-v1.0.0-x86_64-unknown-linux-gnu/trgt && \
        mv trgt-v1.0.0-x86_64-unknown-linux-gnu/trgt /usr/bin/ && \
        rm -r trgt-v1.0.0-x86_64-unknown-linux-gnu/

RUN wget "https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2" && \
        bunzip2 bcftools-1.20.tar.bz2 && \
        tar -xvf bcftools-1.20.tar && \
        cd bcftools-1.20/ && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -r bcftools-1.20/ && \
        rm bcftools-1.20.tar

RUN wget "https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2" && \
        bunzip2 htslib-1.20.tar.bz2 && \
        tar -xvf htslib-1.20.tar && \
        cd htslib-1.20/ && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -r htslib-1.20/ && \
        rm htslib-1.20.tar

RUN wget "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz" && \
        gunzip adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz

RUN wget "https://zuchnerlab.s3.amazonaws.com/RepeatExpansions/TRGT/T2T/CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz" && \
        gunzip CHM13.adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0_under1kb.bed.gz

