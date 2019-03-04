FROM weisburd/base-image-for-str-tools:v1

MAINTAINER Ben Weisburd

ENV VERSION=1.4

# install GangSTR from https://github.com/gymreklab/GangSTR
ENV HTSLIB_CFLAGS="-I./htslib-1.9"
ENV HTSLIB_LIBS=" -L$./htslib-1.9 -lhts"
ENV NLOPT_CFLAGS="-I./nlopt-2.4.2"
ENV NLOPT_LIBS=" -L$./nlopt-2.4.2 -lnlopt"
ENV CPPFLAGS="-std=c++98"

ENV LD_LIBRARY_PATH=/usr/local/lib/

RUN wget https://github.com/gymreklab/GangSTR/releases/download/v${VERSION}/GangSTR-${VERSION}.tar.gz \
    && tar xzf GangSTR-${VERSION}.tar.gz \
    && rm GangSTR-${VERSION}.tar.gz

WORKDIR /GangSTR-${VERSION}

RUN ./install-gangstr.sh

# remove ENTRYPOINT because cromwell doesn't work when it's defined
#ENTRYPOINT [ "GangSTR" ]
