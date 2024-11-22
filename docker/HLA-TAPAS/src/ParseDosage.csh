#!/bin/csh -f

###############################################################
#
# Convert .gprobs to .dos file
# Usage: ./ParseDosage.csh file.gprobs > file.dos
#
###############################################################

set INPUT=$1
awk '{printf("%s\t%s\t%s",$1,$2,$3); for (f = 4; f <= NF; f+=3){printf("\t%.3f",$f*2+$(f+1));} printf("\n");}' $INPUT
