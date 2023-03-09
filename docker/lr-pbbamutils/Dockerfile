FROM continuumio/miniconda3

RUN conda update -n base conda && \
    conda install -n base conda-libmamba-solver && \
    conda config --set solver libmamba

RUN conda install -c conda-forge libgcc-ng=12.2.0 && \
    conda install -c bioconda \
                  pbmm2=1.10.0 \
                  pbbam=2.1.0 \
                  samtools=1.16.1 && \
    conda clean -a

