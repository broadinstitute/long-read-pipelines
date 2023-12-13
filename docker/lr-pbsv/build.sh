#!/bin/bash
set -euo pipefail

conda env create -f environment.yml
conda clean -a
echo "conda activate lr-pbsv" >> ~/.bashrc
