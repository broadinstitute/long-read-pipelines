FROM us.gcr.io/broad-dsp-lrma/lr-papermill-base:2.5.0 AS build

# Install R packages using the R package manager needed by notebook:
RUN R -e "install.packages(c('ggbeeswarm', 'beeswarm', 'dplyr', 'ggplot2', 'ggtext'), repos='https://cran.r-project.org')"

COPY plot_single_sample_stats.ipynb /
