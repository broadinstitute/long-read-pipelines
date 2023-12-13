FROM continuumio/miniconda3
COPY --from=us.gcr.io/broad-dsp-lrma/lrm-gcloud:0.1.1 /install /install
RUN /install/install.sh

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -ay
ENV PATH=/opt/conda/envs/lr-sniffles2/bin/:/root/google-cloud-sdk/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/lr-sniffles2/lib/
RUN echo "source activate lr-sniffles2" >> ~/.bashrc
