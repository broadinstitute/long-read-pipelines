FROM rstudio/r-base:4.4.1-focal

MAINTAINER Yulia Mostovoy

RUN apt-get update -y && sudo apt-get upgrade -y \
    && apt install -y build-essential checkinstall zlib1g-dev libssl-dev libxml2-dev \
    && apt-get clean

RUN Rscript -e "install.packages('optparse', repos='http://cran.rstudio.com/')" \
    && Rscript -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" \
    && Rscript -e "install.packages('remotes', repos='http://cran.rstudio.com/')" \
    && Rscript -e "BiocManager::install(c('BiocParallel', 'edgeR', 'bsseq', 'DMRichR'))"

COPY ./callDMRs.R /
