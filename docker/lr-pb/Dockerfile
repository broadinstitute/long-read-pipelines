FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-pb/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt-get --allow-releaseinfo-change update
RUN apt install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# copy other resources
COPY detect_run_info.py /usr/local/bin/
COPY merge_ccs_reports.py /usr/local/bin/
COPY shard_bam.py /usr/local/bin/
COPY extract_uncorrected_reads.py /usr/local/bin/
COPY collect_zmw_subread_stats.py /usr/local/bin/
COPY collect_polymerase_read_lengths.py /usr/local/bin/
COPY add_pac_bio_rg_to_bam.py /usr/local/bin/
COPY compute_pbi_stats.py /usr/local/bin/
COPY reset_clr_bam_bq.py /usr/local/bin/
COPY compute_read_length_hist.py /usr/local/bin/

# activate conda environment
RUN echo "source activate lr-pb" > ~/.bashrc
