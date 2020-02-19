FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev vim wget unzip \
    && rm -rf /var/lib/apt/lists/*

# install htslib
RUN git clone https://github.com/samtools/htslib.git \
    && cd htslib \
    && autoheader \
    && autoconf \
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

ADD tool.py /lrma/tool.py
ADD tool_rle.py /lrma/tool_rle.py

RUN echo "source activate 10x_tool" > ~/.bashrc
RUN pip3 install pysam biopython
