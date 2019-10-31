FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create python=3.6 -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-longshot/bin/:${PATH}
RUN echo "source activate lr-longshot" > ~/.bashrc
RUN activate lr-longshot
