FROM gpuci/clara-genomics-base:cuda10.1-ubuntu18.04-gcc7-py3.5

RUN apt-get -qqy update && \
    apt-get -qqy install cmake git && \
    cd /tmp && \
    git clone --recursive https://github.com/lbcb-sci/racon.git racon && \
    cd racon && mkdir build && cd build && \
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -Dracon_build_tests=ON \
        -Dracon_enable_cuda=ON \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        .. && \
    make && make install
