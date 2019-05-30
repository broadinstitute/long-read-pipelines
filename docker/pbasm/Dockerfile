FROM kgarimella/pbtools:latest

MAINTAINER Kiran Garimella

ENV PYTHON_VERSION="3.7.1"

RUN git clone https://github.com/lh3/miniasm \
        && cd miniasm \
        && make

# Install Racon (requiring gcc 4.8+ or clang 3.4+, and cmake 3.2+)
RUN git clone --recursive https://github.com/isovic/racon.git racon \
        && cd racon \
        && mkdir build \
        && cd build \
        && cmake -DCMAKE_BUILD_TYPE=Release .. \
        && make

ENV PATH=/miniasm/:/racon/build/bin/:${PATH}
