#!/usr/bin/env bash

input_vcf=$1
output_vcf=$(basename ${input_vcf} | sed s/.vcf/.normalized.vcf/)
fasta=${PROJECT1_BASE}/ref/GRCh37/hg19.fa

echo ${input_vcf} ' > '  ${output_vcf}

gzcat ${input_vcf} | vt decompose -s - | vt normalize -r ${fasta} - | bgzip > ${output_vcf}

${PROJECT1_BASE}/scripts/count_vcf.sh ${output_vcf}
