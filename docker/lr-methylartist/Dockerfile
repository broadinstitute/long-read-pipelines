FROM continuumio/miniconda3
RUN apt-get update
RUN apt-get install -y --no-install-recommends git ssh build-essential make cmake zlib1g-dev vim wget unzip \
    && rm -rf /var/lib/apt/lists/*
    
RUN pip install methylartist