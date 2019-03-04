#!/usr/bin/env bash

input_bed=$1
output_bed=$(basename ${input_bed} | sed s/.bed/.filtered_to_high_confidence_no_polyA.bed/)
bed_file=${PROJECT1_BASE}/data/GRCh37/syndip/CHM-eval.kit/full.37d5.bed.gz

echo ${input_bed} ' / ' ${bed_file} ' > '  ${output_bed}

bedtools intersect -a ${input_bed} -b ${bed_file} -u -f 1  > ${output_bed}

