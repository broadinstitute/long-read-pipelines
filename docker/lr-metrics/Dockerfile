FROM dceoy/r-tidyverse

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install some utilities
RUN apt-get update
RUN apt-get install wget datamash

# install conda packages
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
        && bash /miniconda.sh -b -p /miniconda
ENV PATH=/miniconda/bin/:/miniconda/envs/lr-metrics/bin/:/root/google-cloud-sdk/bin/:${PATH}
RUN conda env create -f /environment.yml && conda clean -a

# install more R packages
RUN set -e \
    && clir update \
    && clir install --devt=cran vroom \
    && clir validate vroom

# remove the entry point that we inherited from the r-tidyverse docker
ENTRYPOINT []
WORKDIR /
CMD ["/bin/bash"]
