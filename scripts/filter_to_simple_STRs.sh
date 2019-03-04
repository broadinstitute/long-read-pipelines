#!/usr/bin/env bash

input_vcf=$1
output_vcf=$(basename ${input_vcf} | sed s/.vcf/.simple_STRs_only.vcf/)

echo filter to simple STRs: ${input_vcf} ' > '  ${output_vcf}

python ${PROJECT1_BASE}/scripts/simple_STR_filter.py $1 | bgzip > ${output_vcf}

${PROJECT1_BASE}/scripts/count_vcf.sh  ${output_vcf}