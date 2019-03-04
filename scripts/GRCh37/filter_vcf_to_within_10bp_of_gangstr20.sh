#!/usr/bin/env bash

input_vcf=$1
output_vcf=$(basename ${input_vcf} | sed s/.vcf/.filtered_to_within_10bp_of_gangstr20.vcf/)
bed_file=${PROJECT1_BASE}/ref/GRCh37/repeats/computed/gangstr20_confidence_region_intersection.padded_10bp.bed

echo ${input_vcf} ' / ' ${bed_file} ' > '  ${output_vcf}

echo Before: 
${PROJECT1_BASE}/scripts/count_vcf.sh ${input_vcf}

bedtools intersect -a ${input_vcf} -b ${bed_file} -u -header | bgzip > ${output_vcf}.gz

echo After:
${PROJECT1_BASE}/scripts/count_vcf.sh ${output_vcf}.gz