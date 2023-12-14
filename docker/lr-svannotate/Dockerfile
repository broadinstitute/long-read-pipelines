FROM broadinstitute/gatk:4.5.0.0

MAINTAINER Yulia Mostovoy

RUN apt-get -y update \
	&& apt-get -y install gcc cmake protobuf-compiler libcurl4-openssl-dev libbz2-dev liblzma-dev bzip2 \
	&& apt-get clean

# download my gatk-sv with custom standardize modules for long-read formats
RUN git clone https://github.com/yuliamostovoy/gatk-sv.git \
	&& cd gatk-sv \
	&& git checkout 4c4cfbd18b97e9bc235b15408118afe596015842 \
	&& git reset --hard \
    && pip install Cython \
    && pip install -e src/svtk/

# add additional processing script
COPY format_svtk_vcf_for_gatk.py .
