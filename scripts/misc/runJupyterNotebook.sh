#!/usr/bin/env bash

# This script is designed to be run on a machine on which you want to start up a jupyter lab server.
#
# Author: Jonn Smith

d=$1
if [[ ${#d} -eq 0 ]] ; then
	d="${PWD}"
fi

docker run -it -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v ${d}:/home/jovyan/work us.gcr.io/broad-dsp-lrma/lr-jupyter:latest 

