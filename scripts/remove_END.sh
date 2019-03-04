export F=$1
#export F=/Users/weisburd/project1/results2/gangstr/gangstr_syndip_genomes_CHM1_CHM13_2.standard_format.vcf.gz; 

gzcat $F | sed 's/END\=[0-9]*;//' | bgzip > ${F}.temp.vcf.gz && mv ${F}.temp.vcf.gz ${F}
