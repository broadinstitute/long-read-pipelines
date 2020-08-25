FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# install some utilities
RUN apt-get update -m
RUN apt-get -y install bsdmainutils git-lfs wget datamash autoconf make cmake gcc g++ zlib1g-dev
RUN wget http://ftp.de.debian.org/debian/pool/contrib/m/msttcorefonts/ttf-mscorefonts-installer_3.7_all.deb \
    && apt-get -y install /ttf-mscorefonts-installer_3.7_all.deb

# install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
    && bash /miniconda.sh -b -p /miniconda
ENV PATH=/miniconda/bin/:/miniconda/envs/lr-metrics/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# install gatk
RUN git clone https://github.com/broadinstitute/gatk.git -b kvg_pbeap \
    && cd gatk \
    && ./gradlew \
    && cp /gatk/build/libs/gatk-package-4.1.0.0-54-g47086b6-SNAPSHOT-local.jar /usr/local/bin/gatk.jar \
    && cd .. \
    && rm -rf /gatk

# install picard
RUN wget -O /usr/local/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/2.22.1/picard.jar

# install various metric and visualization scripts
COPY lima_report_detail.R /
COPY lima_report_summary.R /
COPY prl_plots.R /
COPY utils.R /
COPY yield.R /
COPY Sequel_96_barcodes_v2.fasta /
