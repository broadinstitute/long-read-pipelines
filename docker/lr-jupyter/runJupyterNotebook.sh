#!/usr/bin/env bash

d=$1
if [[ ${#d} -eq 0 ]] ; then
	d="${PWD}"
fi

#docker run -it -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v ${d}:/home/jovyan/work jonnsmith/jupyternotebook:latest
docker run -it -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v ${d}:/home/jovyan/work us.gcr.io/broad-dsp-lrma/lr-jupyter:latest 

