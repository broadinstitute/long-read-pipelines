PWD=`pwd`

.PHONY : all

pbeap-base/sha256.txt: pbeap-base/Dockerfile
	printf "\n============\nPBEAP Base Image\n" \
	&& cd $(PWD)/pbeap-base && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt

pbtools/sha256.txt: pbeap-base pbtools/Dockerfile
	printf "\n============\nPB Tools Image\n" \
	&& cd $(PWD)/pbtools && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt

all: pbtools
