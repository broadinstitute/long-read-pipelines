# Start with a good base python3 image:
FROM ubuntu:20.04
MAINTAINER Jonn Smith <jonn@broadinstitute.org>

# Make sure we don't need to interact with any package installations:
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /
WORKDIR /

########################################################################################################################
# DEPENDENCIES

# Install some dependencies for gsutil:
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git-lfs parallel

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3 python3-pip python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# Install gsutil:
RUN curl https://sdk.cloud.google.com | bash

# Dependencies for samtools:
RUN apt-get install -y bzip2 curl gnupg2 libc-dev ncurses-dev libcurl4-openssl-dev libssl-dev libbz2-dev liblzma-dev

# Google cloud support for samtools:
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
		&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
		&& apt-get update -y \
		&& apt-get install google-cloud-sdk -y

# Additional Dependencies:
RUN apt install -y curl wget datamash pkg-config zip unzip default-jre python

########################################################################################################################
# SOFTWARE:

RUN mkdir -p /usr/local/bin /usr/local/lib /usr/local/etc

# BWA:
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwakit-0.7.15_x64-linux.tar.bz2 && \
    tar -xf bwakit-0.7.15_x64-linux.tar.bz2 && \
    cd bwa.kit && \
    mv bwa fermi2 fermi2.pl htsbox k8 ropebwt2 run-HLA run-bwamem run-gen-ref samblaster samtools seqtk trimadap typeHLA.sh typeHLA-selctg.js typeHLA.js bwa-postalt.js /usr/local/bin && \
    mkdir /usr/local/etc/bwa && \
    mv README.md doc resource-GRCh38 resource-human-HLA /usr/local/etc/bwa && \
    cd .. && \
    rm -rf bwakit-0.7.15_x64-linux.tar.bz2 bwa.kit

# TABIX:
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
    tar -xf htslib-1.18.tar.bz2 && \
    cd htslib-1.18 && \
    ./configure && make && make install && \
    cd .. \
    rm -rf htslib-1.18 htslib-1.18.tar.bz2

# VCF Tools:
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure && make && make install && \
    cd .. \
    rm -rf vcftools-0.1.16 vcftools-0.1.16.tar.gz

# bcftools:
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -xf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && \
    ./configure && make && make install && \
    cd .. \
    rm -rf bcftools-1.9 bcftools-1.9.tar.bz2

# Samtools:
# Get samtools source:
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar -xjf samtools-1.11.tar.bz2 && \
	cd samtools-1.11 && \
	./configure && make install && \
    cd .. && \
    rm -rf samtools-1.11 samtools-1.11.tar.bz2

# gatk/4.2.2.0:
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip && \
    unzip gatk-4.2.2.0.zip && \
    cd gatk-4.2.2.0 && \
    mv gatk gatk-completion.sh gatk-package-4.2.2.0-local.jar gatk-package-4.2.2.0-spark.jar /usr/local/bin && \
    mkdir -p /usr/local/etc/gatk && \
    mv GATKConfig.EXAMPLE.properties README.md gatkPythonPackageArchive.zip gatkcondaenv.yml gatkdoc scripts /usr/local/etc/gatk && \
    cd .. \
    rm -rf gatk-4.2.2.0.zip gatk-4.2.2.0

# zlib/1.2.11:
RUN wget https://www.zlib.net/fossils/zlib-1.2.11.tar.gz && \
    tar -xf zlib-1.2.11.tar.gz && \
    cd zlib-1.2.11 && \
    ./configure && make && make install && \
    cd .. \
    rm -rf zlib-1.2.11 zlib-1.2.11.tar.gz

# plink/1.90:
RUN mkdir plink && \
    cd plink && \
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip && \
    unzip plink_linux_x86_64_20231018.zip && \
    rm plink_linux_x86_64_20231018.zip && \
    mv plink prettify /usr/local/bin && \
    mkdir /usr/local/etc/plink && \
    mv toy.ped toy.map /usr/local/etc/plink/ && \
    cd .. && \
    rm -rf plink

# trimmomatic/0.36:
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
    mkdir -p /opt && \
    unzip Trimmomatic-0.36.zip && \
    mv Trimmomatic-0.36 /opt/ && \
    echo "#!/usr/bin/env bash\n\njava -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE $@\n\n" > /usr/local/bin/TrimmomaticPE && \
    echo "#!/usr/bin/env bash\n\njava -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar SE $@\n\n" > /usr/local/bin/TrimmomaticSE && \
    chmod +x /usr/local/bin/Trimmomatic*E && \
    rm Trimmomatic-0.36.zip

########################################################################################################################
########################################################################################################################
########################################################################################################################

# Other utilities:
RUN apt-get clean
RUN apt install -y vim emacs nano
RUN apt-get clean

########################################################################################################################
########################################################################################################################
########################################################################################################################

# Might not need these:

# gsl/2.7.1:

# RAiSD/2.8:

# R/4.1.0:

# sratoolkit/2.8.2-1:
# BROKEN
#RUN wget https://github.com/ncbi/sra-tools/archive/refs/tags/2.8.2-1.tar.gz && \
#    tar -xf 2.8.2-1.tar.gz && \
#    cd sra-tools-2.8.2-1 && \


