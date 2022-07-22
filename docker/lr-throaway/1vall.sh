#!/bin/bash

set -eux

aligned_bam=$1
fingerprint_vcf=$2
haplotype_map=$3
prefix=$4

#################################################
GREPCMD="grep"
if [[ "${fingerprint_vcf}" =~ \.gz$ ]]; then
    GREPCMD="zgrep"
fi

vcf_sample_name=$("${GREPCMD}" "^#CHROM" "${fingerprint_vcf}" | awk '{print $10}')

#################################################
set +e
gatk CheckFingerprint \
    --INPUT "${aligned_bam}" \
    --GENOTYPES "${fingerprint_vcf}" \
    --EXPECTED_SAMPLE_ALIAS "${vcf_sample_name}" \
    --HAPLOTYPE_MAP "${haplotype_map}" \
    --OUTPUT "${prefix}"
# non-zero status indicates something wrong, touch up a dummy map and exit
if [[ $? -ne 0 ]]; then
    echo "CheckFingerprint failed. Patchy work then exit 0."
    echo -e "LOD_EXPECTED_SAMPLE\t-10000.0" > "${prefix}.metrics_map.txt"
    exit 0;
fi
set -e
#################################################
grep -v '^#' "${prefix}".fingerprinting_summary_metrics | \
    grep -A1 READ_GROUP | \
    awk '
        {
            for (i=1; i<=NF; i++)  {
                a[NR,i] = $i
            }
        }
        NF>p { p = NF }
        END {
            for(j=1; j<=p; j++) {
                str=a[1,j]
                for(i=2; i<=NR; i++){
                    str=str" "a[i,j];
                }
                print str
            }
        }' | \
    sed 's/ /\t/' \
    > "${prefix}.metrics_map.txt"

mv "${prefix}.fingerprinting_summary_metrics" \
    "${prefix}.fingerprinting_summary_metrics.txt"
mv "${prefix}.fingerprinting_detail_metrics" \
    "${prefix}.fingerprinting_detail_metrics.txt"
