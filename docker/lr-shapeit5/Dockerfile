FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/lr-shapeit5/bin/:${PATH}

# install packages
RUN apt-get update --allow-releaseinfo-change
RUN apt-get install -y --no-install-recommends git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev curl git-lfs time gnupg autoconf bzip2 bcftools \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/*

# install SHAPEIT5
RUN wget -P /bin/ https://github.com/odelaneau/shapeit5/releases/download/v1.0.0-beta/SHAPEIT5_phase_common_static
RUN wget -P /bin/ https://github.com/odelaneau/shapeit5/releases/download/v1.0.0-beta/SHAPEIT5_phase_rare_static
RUN wget -P /bin/ https://github.com/odelaneau/shapeit5/releases/download/v1.0.0-beta/SHAPEIT5_switch_static
RUN wget -P /bin/ https://github.com/odelaneau/shapeit5/releases/download/v1.0.0-beta/SHAPEIT5_ligate_static
RUN wget -P /bin/ https://github.com/odelaneau/shapeit5/releases/download/v1.0.0-beta/SHAPEIT5_convert_static

RUN chmod +x /bin/SHAPEIT5*
