input_vcf=$1

#echo input: ${input_vcf} 

python ${PROJECT1_BASE}/scripts/vcf_variant_size_distribution.py $1 

#${PROJECT1_BASE}/scripts/count_vcf.sh  ${output_vcf}