FROM r-base:4.3.0

RUN Rscript -e "update.packages(contriburl=contrib.url(repos=c('http://lib.stat.cmu.edu/R/CRAN/', 'https://cran.rstudio.com')), clean=TRUE, quiet=TRUE, ask=FALSE)"
COPY install_R_packages.R /opt/
ARG R_PKGS="stringr ggplot2 reshape2"
RUN chmod +x /opt/install_R_packages.R && \
    /opt/install_R_packages.R $R_PKGS

COPY plot.resources.R /opt/
RUN chmod +x /opt/plot.resources.R
