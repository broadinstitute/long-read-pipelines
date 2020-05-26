FROM gpuci/clara-genomics-base:cuda10.1-ubuntu18.04-gcc7-py3.5

ARG RACON_REPO="https://github.com/lbcb-sci/racon.git"
ARG RACON_BRANCH="master"
RUN apt-get -qqy update && \
    apt-get -qqy install cmake git && \
    cd /tmp && \
    git clone --recursive $RACON_REPO racon && \
    cd racon && git checkout $RACON_BRANCH && \
    mkdir build && cd build && \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -Dracon_build_tests=ON \
        -Dracon_enable_cuda=ON \
        -Dracon_build_wrapper=ON \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        .. && \
    make && make install

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -  && \
    mv minimap2-2.17_x64-linux/minimap2 /usr/local/bin