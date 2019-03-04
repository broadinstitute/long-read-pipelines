input_vcf=$1
output_vcf=$(basename ${input_vcf} | sed s/.vcf/.deletions.vcf/)

echo select deletions in ${input_vcf} ' > '  ${output_vcf}

cat <(bcftools view -h ${input_vcf}) <(bcftools view -v indels ${input_vcf} | awk 'length($4) > length($5) { print }') | bgzip > ${output_vcf}

${PROJECT1_BASE}/scripts/count_vcf.sh ${output_vcf}
