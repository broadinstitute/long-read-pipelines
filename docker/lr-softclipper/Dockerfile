FROM continuumio/miniconda3

RUN conda install python=3.6.9
RUN conda create -n venv python=3.6.9
RUN activate venv

RUN conda install -c bioconda samtools=1.9
RUN conda install -c bioconda minimap2

COPY soft_clipper.py .
