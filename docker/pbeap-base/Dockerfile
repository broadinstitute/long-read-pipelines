# This base image contains dependencies used by the PacBio EAP processing pipeline.
# It uses bitnami/minideb:stretch as its base image to work with cromwell.
# Most of this is lifted from Ben Weisburd's STR analysis repository (https://github.com/bw2/str-callers).

FROM bitnami/minideb:stretch

MAINTAINER Kiran Garimella

ENV PYTHON_VERSION="3.7.1"
ENV GATK_VERSION="4.1.0.0"
ENV PICARD_VERSION="2.18.27"
ENV SAMTOOLS_VERSION="1.9"
ENV VCF_TOOLS_VERSION="0.1.16"
ENV BED_TOOLS_VERSION="2.27.1"

RUN install_packages ca-certificates wget pkg-config gcc g++ make bzip2 zlib1g-dev dpkg-dev build-essential libcurl4-openssl-dev libbz2-dev liblzma-dev

# python2 utils
RUN install_packages python-dev python-setuptools python-pip \
	&& pip install --upgrade wheel

# install python3.8
RUN apt-get update \
	&& apt-get dist-upgrade -y \
	&& apt-get install -y python-smbus libncursesw5-dev libgdbm-dev libc6-dev zlib1g-dev libsqlite3-dev libssl-dev openssl libffi-dev

RUN wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
	&& tar xvf Python-${PYTHON_VERSION}.tar.xz \
	&& rm Python-${PYTHON_VERSION}.tar.xz \
	&& cd Python-${PYTHON_VERSION} \
	&& ./configure --enable-optimizations \
	&& make install

# install java-8
RUN apt install -y default-jdk

# install picard
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -O /picard.jar

# install GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
	&& apt install -y unzip \
	&& unzip gatk-${GATK_VERSION}.zip \
	&& rm gatk-${GATK_VERSION}.zip \
	&& ln -s gatk-${GATK_VERSION}/gatk-package-${GATK_VERSION}-local.jar /gatk.jar

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install

#RUN apt install -y autoconf unzip \
#	&& wget https://github.com/samtools/bcftools/archive/develop.zip -O bcftools-develop.zip \
#	&& unzip bcftools-develop.zip \
#	&& rm bcftools-develop.zip \
#	&& cd bcftools-develop \
#	&& autoconf \
#	&& autoheader \
#	&& ./configure \
#	&& make \
#	&& make install

# install vcf-tools
RUN wget https://github.com/vcftools/vcftools/releases/download/v${VCF_TOOLS_VERSION}/vcftools-${VCF_TOOLS_VERSION}.tar.gz \
	&& tar xzf vcftools-${VCF_TOOLS_VERSION}.tar.gz \
	&& rm vcftools-${VCF_TOOLS_VERSION}.tar.gz \
	&& cd vcftools-${VCF_TOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install

# install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BED_TOOLS_VERSION}/bedtools-${BED_TOOLS_VERSION}.tar.gz \
	&& tar xzf bedtools-${BED_TOOLS_VERSION}.tar.gz \
	&& rm bedtools-${BED_TOOLS_VERSION}.tar.gz \
	&& cd bedtools2 \
	&& make \
	&& make install

# install verticalize (https://github.com/lindenb/verticalize)
RUN wget https://github.com/lindenb/verticalize/archive/v1.0.1.tar.gz -O verticalize-1.0.1.tar.gz \
	&& tar xzf verticalize-1.0.1.tar.gz \
	&& cd verticalize-1.0.1 \
	&& make \
	&& mv verticalize /usr/local/bin \
	&& cd .. \
	&& rm -rf verticalize-1.0.1*

# igvtools
RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.4/igvtools_2.4.17.zip \
	&& unzip igvtools_2.4.17

# nice-to-have linux utilities
RUN install_packages emacs vim less curl htop xterm git

# python3 packages
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install --upgrade tqdm pyvcf pyfaidx cyvcf2

# miniconda installation
RUN curl -LO https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/root/google-cloud-sdk/bin/:/miniconda/bin:${PATH}
RUN conda update -y conda

# configure conda
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels r
RUN conda update --all

# install conda packages
RUN conda install -c bioconda pbccs
RUN conda install -c bioconda pbmm2
RUN conda install -c bioconda pbsv
RUN conda install -c bioconda nanoplot
RUN conda install -c bioconda fastqc
RUN conda install -c bioconda bcftools
RUN conda install -c bioconda vcftools
RUN conda install -c r r

# script for splitting .bed file into multiple .bed files
COPY scatter_intervals.py /
COPY scatter_intervals_by_covered_region_size.py /
COPY compute_print_read_intervals.py /
COPY generate_repeat_specs.py /
COPY standardize_output_vcf.py /
COPY update_vcf_id_field.py /

COPY bashrc /root/.bashrc

# install gsutil
RUN curl https://sdk.cloud.google.com | bash
#ENV PATH=/root/google-cloud-sdk/bin/:${PATH}
#RUN gcloud init
