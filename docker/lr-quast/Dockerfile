FROM continuumio/miniconda3

RUN conda install python=3.6.9
RUN conda create -n lr-quast python=3.6.9
RUN activate lr-quast

RUN conda install matplotlib
RUN conda install -c bioconda joblib
RUN conda install -c bioconda quast