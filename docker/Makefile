PWD=`pwd`

.PHONY : pbeap-base

pbeap-base:
	printf "\n============\nPBEAP Base Image\n" \
	&& cd $(PWD)/pbeap-base && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt

#gangstr: base-image
#	printf "\n============\n GangSTR \n" \
#	&& cd $(PWD)/gangstr && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
#	&&   sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/gangstr/*.wdl \
#	&& cat sha256.txt

all: pbeap-base
