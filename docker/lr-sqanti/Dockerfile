FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

RUN apt-get -y update \
	&& apt-get -y install git make gcc g++ zlib1g-dev tree wget \
	&& apt-get clean

# Install SQANTI3 and prepare Python environment
RUN git clone https://github.com/ConesaLab/SQANTI3.git \
    && cd SQANTI3 \
    && conda env create -f SQANTI3.conda_env.yml \
    && /bin/bash -c "source activate SQANTI3.env"

# Install gtfToGenePred
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P /SQANTI3/utilities/ \
    && chmod +x SQANTI3/utilities/gtfToGenePred 

# Install cDNA_Cupcake
RUN git clone https://github.com/Magdoll/cDNA_Cupcake.git \
    && conda install numpy cython \
    && cd cDNA_Cupcake \
    && python setup.py build \
    && python setup.py install

# Configure environment
ENV PATH=${PATH}:/SQANTI3/utilities/
ENV PYTHONPATH=${PYTHONPATH}:/cDNA_Cupcake/sequence/
ENV PYTHONPATH=${PYTHONPATH}:/cDNA_Cupcake/
RUN echo "source activate SQANTI3.env" > ~/.bashrc
