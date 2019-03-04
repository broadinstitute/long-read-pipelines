FROM weisburd/base-image-for-str-tools:v1

MAINTAINER Ben Weisburd

RUN wget http://tandem.bu.edu/trf/downloads/trf409.linux64 \
	&& chmod 777 /trf409.linux64

RUN python3 -m pip install --upgrade peewee

COPY convert_dat_output_to_bed.py /
COPY dat_utils.py /
COPY bed_utils.py /
COPY filter_vcf_to_STR_alleles.py /
COPY export_dat_or_bed_to_database.py /