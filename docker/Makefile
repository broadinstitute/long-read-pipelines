MODULES = lr-align lr-asm lr-bri lr-c3poa lr-cloud-downloader lr-finalize lr-gatk lr-guppy lr-longshot lr-metrics lr-ont lr-pb lr-peregrine lr-sv lr-svim lr-utils racon sr-utils

all:
	for dir in $(MODULES); do \
		echo begin $$dir; \
		(cd $$dir; ${MAKE}); \
		echo end $$dir; \
		echo ==========; \
	done
