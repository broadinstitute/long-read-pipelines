#!/usr/bin/env bash

input_vcf=$1
output_vcf=$(basename ${input_vcf} | sed s/.vcf/.filtered_to_high_confidence_no_polyA.vcf/)
bed_file=${PROJECT1_BASE}/data/GRCh38/syndip/CHM-eval.kit/full.38.bed.gz

echo ${input_vcf} ' / ' ${bed_file} ' > '  ${output_vcf}

bedtools intersect -a ${input_vcf} -b ${bed_file} -u -f 1 -header | bgzip > ${output_vcf}

${PROJECT1_BASE}/scripts/count_vcf.sh ${output_vcf}