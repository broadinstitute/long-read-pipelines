FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev vim wget unzip \
    && rm -rf /var/lib/apt/lists/*

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 \
    && tar -jxf htslib-1.11.tar.bz2 \
    && cd htslib-1.11 \
    && ./configure \
    && make \
    && make install

# install bri
RUN git clone https://github.com/jts/bri.git \
    && cd bri \
    && make \
    && cp bri /usr/local/bin/
ENV LD_LIBRARY_PATH=/usr/local/lib/

WORKDIR /lrma

# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
    && bash /miniconda.sh -b -p /miniconda

# Install samtools 1.11:
RUN apt-get update
RUN apt-get install -y bzip2 curl gnupg2
RUN apt-get install -y libc-dev ncurses-dev

#### Specific for google cloud support
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
		&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - \
		&& apt-get update -y \
		&& apt-get install google-cloud-sdk -y

RUN apt-get install -y libcurl4-openssl-dev libssl-dev

# Get samtools source:
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
		&& tar -xjf samtools-1.11.tar.bz2 \
		&& cd samtools-1.11 \
		&& ./configure \
		&& make install
RUN rm -rf /samtools-1.11 /samtools-1.11.tar.bz2

# install conda packages
ADD ./environment.yml /lrma/environment.yml
RUN pip install --upgrade pip
RUN conda env create -f environment.yml

# install annotation software
ADD ./ssw_lib.py /lrma/ssw/ssw_lib.py
ADD ./__init__.py /lrma/ssw/__init__.py
RUN wget -O ssw-zip.zip https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/master.zip
RUN unzip ssw-zip.zip
RUN cd Complete-Striped-Smith-Waterman-Library-master/src && make libssw.so && cd ../..
RUN mv Complete-Striped-Smith-Waterman-Library-master/src/libssw.so /lrma/ssw/libssw.so

RUN wget -O starcode-zip.zip https://github.com/gui11aume/starcode/archive/master.zip
RUN unzip starcode-zip.zip
RUN cd starcode-master && make && cd ..

ADD 3M-february-2018.txt.gz /lrma/3M-february-2018.txt.gz
RUN gunzip -c /lrma/3M-february-2018.txt.gz > /lrma/3M-february-2018.txt
ADD adapter_sequence.fasta /lrma/adapter_sequence.fasta
ADD adapter_sequence.fasta.amb /lrma/adapter_sequence.fasta.amb
ADD adapter_sequence.fasta.ann /lrma/adapter_sequence.fasta.ann
ADD adapter_sequence.fasta.bwt /lrma/adapter_sequence.fasta.bwt
ADD adapter_sequence.fasta.pac /lrma/adapter_sequence.fasta.pac
ADD adapter_sequence.fasta.sa /lrma/adapter_sequence.fasta.sa
ADD reverse_adapter_sequence.fasta /lrma/reverse_adapter_sequence.fasta
ADD reverse_adapter_sequence.fasta.amb /lrma/reverse_adapter_sequence.fasta.amb
ADD reverse_adapter_sequence.fasta.ann /lrma/reverse_adapter_sequence.fasta.ann
ADD reverse_adapter_sequence.fasta.bwt /lrma/reverse_adapter_sequence.fasta.bwt
ADD reverse_adapter_sequence.fasta.pac /lrma/reverse_adapter_sequence.fasta.pac
ADD reverse_adapter_sequence.fasta.sa /lrma/reverse_adapter_sequence.fasta.sa

RUN pip3 install pysam biopython
RUN pip3 install tqdm

ADD tool.py /lrma/tool.py
ADD tool_rle.py /lrma/tool_rle.py
ADD tool_starcode_seeded.py /lrma/tool_starcode_seeded.py
ADD extract_ilmn_bc_conf_scores.py /lrma/extract_ilmn_bc_conf_scores.py

ADD restore_annotations_to_aligned_bam.py /lrma/restore_annotations_to_aligned_bam.py
ADD extract_cbc_and_umi_from_annotated_read.py /lrma/extract_cbc_and_umi_from_annotated_read.py 
ADD copy_contig_name_to_tag.py /lrma/copy_contig_name_to_tag.py 
ADD tag_mas_sirv_umi_positions.py /lrma/tag_mas_sirv_umi_positions.py 

RUN echo "source activate 10x_tool" > ~/.bashrc

