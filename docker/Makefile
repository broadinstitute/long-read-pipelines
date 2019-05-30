PWD=`pwd`

.PHONY : pbeap-base pbtools pbasm

pbtools: pbeap-base pbtools/Dockerfile
	printf "\n============\nPB Tools Image\n" \
	&& cd $(PWD)/pbtools && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt

pbeap-base: pbeap-base/Dockerfile
	printf "\n============\nPBEAP Base Image\n" \
	&& cd $(PWD)/pbeap-base && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt

pbasm: pbtools pbasm/Dockerfile
	printf "\n============\nPBEAP ASM\n" \
	&& cd $(PWD)/pbasm && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt

all: pbeap-base pbtools pbasm
