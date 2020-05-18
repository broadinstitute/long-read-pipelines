FROM continuumio/miniconda3

RUN apt-get install -y build-essential make gcc libz-dev libncurses5-dev libbz2-dev zlib1g-dev liblzma-dev libcurl4-openssl-dev parallel curl

RUN curl https://sdk.cloud.google.com | bash
ENV PATH="/root/google-cloud-sdk/bin:${PATH}"

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar -jxvf samtools-1.10.tar.bz2
WORKDIR samtools-1.10
RUN make install
WORKDIR /

RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && \
   tar -jxvf minimap2-2.17_x64-linux.tar.bz2  && \
   ln -s /minimap2-2.17_x64-linux/minimap2 /usr/local/bin

RUN git clone --recursive https://github.com/jts/nanopolish.git

WORKDIR nanopolish
RUN make
RUN ln -s /nanopolish/nanopolish /usr/local/bin/

RUN pip install biopython

WORKDIR /home/root