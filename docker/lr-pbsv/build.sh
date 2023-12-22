#!/bin/bash
set -euo pipefail

apt-get update --fix-missing
export DEBIAN_FRONTEND=noninteractive
apt-get -qqy install wget

mkdir -p /miniconda3
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda3/miniconda.sh
bash /miniconda3/miniconda.sh -b -u -p /miniconda3
rm -rf /miniconda3/miniconda.sh
/miniconda3/bin/conda init bash

/miniconda3/bin/conda install -c bioconda pbsv=2.9.0

#this shouldn't be necessary, but I'm getting errors
ln -s /miniconda3/bin/pbsv /usr/local/bin/pbsv

apt-get -qqy purge wget
apt-get -qqy autoremove
apt-get -qqy clean
rm -rf /tmp/* \
       /var/tmp/* \
       /var/cache/apt/* \
       /var/lib/apt/lists/* \
       /usr/share/man/?? \
       /usr/share/man/??_*
