#!/usr/bin/env bash

input_vcf=$1
output_vcf=$(basename ${input_vcf} | sed s/.vcf/.filtered_to_gangstr20.vcf/)
bed_file=${PROJECT1_BASE}/ref/GRCh37/repeats/computed/gangstr20_confidence_region_intersection.bed

echo ${input_vcf} ' / ' ${bed_file} ' > '  ${output_vcf}

bedtools intersect -a ${input_vcf} -b ${bed_file} -u -header | bgzip > ${output_vcf}

${PROJECT1_BASE}/scripts/count_vcf.sh ${output_vcf}