FROM continuumio/miniconda3
RUN apt-get update
RUN apt-get install -y --no-install-recommends git ssh build-essential make cmake zlib1g-dev vim wget unzip \
    && rm -rf /var/lib/apt/lists/*
    
RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.2/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz \
&& tar -xzf pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz \
&& pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores --help

ENV PATH="pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/:$PATH"