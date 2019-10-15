FROM cschin/peregrine:0.1.5.3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-peregrine/bin/:/opt/conda/envs/peregrine/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install some utilities
RUN apt-get install tree

# install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 \
        && tar xjf minimap2-2.17_x64-linux.tar.bz2 \
        && rm minimap2-2.17_x64-linux.tar.bz2 \
        && cd minimap2-2.17_x64-linux \
        && cp minimap2 /usr/local/bin/

# reset entrypoint
ENTRYPOINT []
CMD ["/bin/bash"]
