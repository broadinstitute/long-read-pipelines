PWD=`pwd`

.PHONY : base-image gangstr expansion-hunter expansion-hunter3 gatk-str stretch

base-image:
	printf "\n============\n Base Image \n" \
	&& cd $(PWD)/base-image-for-str-tools && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/other/*.wdl \
	&& cat sha256.txt

gangstr: base-image
	printf "\n============\n GangSTR \n" \
	&& cd $(PWD)/gangstr && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&&   sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/gangstr/*.wdl \
	&& cat sha256.txt

expansion-hunter: base-image
	printf "\n============\n ExpansionHunter \n" \
	&& cd $(PWD)/expansionhunter && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/expansionhunter/*.wdl \
	&& cat sha256.txt \

expansion-hunter3: base-image
	printf "\n============\n ExpansionHunter3 \n" \
	&& cd $(PWD)/expansionhunter3 && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/expansionhunter3/*.wdl \
	&& cat sha256.txt

gatk-str: base-image
	printf "\n============\n Gatk-Str \n" \
	&& cd $(PWD)/gatk-str && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/gatk-str/*.wdl \
	&& cat sha256.txt

stretch: base-image
	printf "\n============\n Stretch \n" \
	&& cd $(PWD)/stretch && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/stretch/*.wdl \
	&& cat sha256.txt

download-pacbio-reads: base-image
	printf "\n============\n Download PacBio Reads \n" \
	&& cd $(PWD)/download-pacbio-reads && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/download-pacbio-reads/*.wdl \
	&& cat sha256.txt

repeat-finder: base-image
	printf "\n============\n Repeat Finder \n" \
	&& cd $(PWD)/repeat-finder && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& sed -i.bak "s/sha256:[^\"]*/"`cat sha256.txt`"/"  ../../wdl/repeat-finder/*.wdl \
	&& cat sha256.txt


all: base-image repeat-finder gangstr expansion-hunter expansion-hunter3 gatk-str stretch download-pacbio-reads 
