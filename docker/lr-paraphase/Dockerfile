FROM continuumio/miniconda3

RUN conda update -n base conda && \
    conda install -n base conda-libmamba-solver && \
    conda config --set solver libmamba

RUN conda install -c conda-forge -c bioconda \
                  paraphase=1.1.3-0 \
                  samtools=1.16.1 \
                  minimap2=2.24
