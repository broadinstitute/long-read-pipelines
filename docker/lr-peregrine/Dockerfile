FROM cschin/peregrine:0.1.5.3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-asm/bin/:/root/google-cloud-sdk/bin/:${PATH}

# reset entrypoint
ENTRYPOINT ["/bin/bash"]
