FROM weisburd/base-image-for-str-tools:v1

MAINTAINER Ben Weisburd

# based on https://github.com/Oshlack/STRetch/wiki/Installing-STRetch
ENV VERSION=0.3.4

# clone STRetch
RUN wget "http://github.com/Oshlack/STRetch/archive/v${VERSION}.tar.gz" \
    && tar xzfv v${VERSION}.tar.gz \
    && rm v${VERSION}.tar.gz

WORKDIR /STRetch-${VERSION}

RUN ln -s /STRetch-${VERSION} /STRetch


# prevent install.sh from downloading the 8Gb reference data zip - this should be downloaded outside the docker image
RUN mkdir -p reference-data && touch reference-data/hg19.PCRfreeWGS_143_STRetch_controls.tsv

RUN ./install.sh

RUN rm -rf reference-data/

# install bwa
RUN wget https://github.com/lh3/bwa/archive/master.zip -O bwa.zip \
	&& unzip bwa.zip \
	&& cd bwa-master \
	&& make \
	&& cp bwa /STRetch/tools/bin/ \
	&& cd .. \
	&& rm -rf bwa*

COPY PCRfreeWGS.controls.tsv /
COPY pipeline_config.groovy ./pipelines/
COPY scripts/*.py ./scripts/
