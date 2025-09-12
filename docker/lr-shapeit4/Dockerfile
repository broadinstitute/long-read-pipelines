FROM continuumio/miniconda3
ARG work_dir=/shapeit4
WORKDIR ${work_dir}


# --------------------------------- Versions -----------------------------------

ARG htslib_version=1.19
ARG samtools_version=1.9
ARG bcftools_version=1.18
# ------------------------------------------------------------------------------

RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        build-essential \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        python-dev-is-python3 \
        bsdmainutils \
        unzip \
        python3-pycurl \
        bedtools \
        libboost-all-dev 
                
                

# HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    ./configure --enable-libcurl --enable-s3 --enable-gcs\
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help





# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && make -s all all-htslib \
    && make install install-htslib \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

ENV LD_LIBRARY_PATH=/usr/local/include/htslib

RUN git clone https://github.com/odelaneau/shapeit4.git \
    &&     cd shapeit4 \
    &&     cat makefile | sed 's|HTSLIB_INC=$(HOME)/Tools/htslib-1.15|HTSLIB_INC=/usr/local/include|' | sed 's|HTSLIB_LIB=$(HOME)/Tools/htslib-1.15/libhts.a|HTSLIB_LIB=/usr/local/lib/libhts.a|' > makefile2 \
    &&     mv makefile makefile.old \
    &&     mv makefile2 makefile  && \
    # sed -i 's|HTSLIB_INC=.*|HTSLIB_INC=/usr/local/include/|' makefile && \
	# sed -i 's|HTSLIB_LIB=.*|HTSLIB_LIB=/usr/local/lib/libhts.a|' makefile && \
	# sed -i 's|BOOST_INC=.*|BOOST_INC=/usr/include/|' makefile && \
	# sed -i 's|BOOST_LIB_IO=.*|BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a|' makefile && \
	# sed -i 's|BOOST_LIB_PO=.*|BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a|' makefile && \
	sed -i 's/DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma -lcurl -lssl -lcrypto/DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma -lcurl -lssl -lcrypto -ldl -lhts -ldeflate/' makefile && \
	make 
    # && make install && make clean  

ENV PATH=/shapeit4/shapeit4/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# # BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --help



