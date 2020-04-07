FROM cschin/peregrine:0.1.6.1

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/peregrine/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install some utilities
RUN apt-get -y install tree

# reset entrypoint
ENTRYPOINT []
CMD ["/bin/bash"]
