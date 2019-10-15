#!/bin/bash

set -euo pipefail

IN_VCF="${1}"
OUT_VCF="${2}"
REF="${3}"
GATK="/gatk/gatk"

echo "Extracting and filtering SNPs from variant call file ${IN_VCF}."
# https://gatkforums.broadinstitute.org/gatk/discussion/11097/cant-use-vqsr-on-non-model-organism-or-small-dataset
# omitting MQ, MQRankSum because reads were filtered at MAPQ60
# omitting InbreedingCoeff because we only have one sample
# omitting FS and SOR because CCS represents information from both strands
# changes 1bp indel QD filter to 5.0

prefix=$(echo "${IN_VCF}" | sed 's/.vcf.gz//')
MULTIALLELIC="${prefix}.multiallelic.raw.vcf"
DECOMPOSED="${prefix}.multiallelic.decomposed.raw.vcf"

DECOMPOSED_SNP="${prefix}.multiallelic.SNP.raw.vcf"
DECOMPOSED_1BP_INDEL="${prefix}.multiallelic.1bpINDEL.raw.vcf"
DECOMPOSED_gt1BP_INDEL="${prefix}.multiallelic.gt1bpINDEL.raw.vcf"

FILTERED_DECOMPOSED_SNP="${prefix}.multiallelic.SNP.filtered.vcf"
FILTERED_DECOMPOSED_1BP_INDEL="${prefix}.multiallelic.1bpINDEL.filtered.vcf"
FILTERED_DECOMPOSED_gt1BP_INDEL="${prefix}.multiallelic.gt1bpINDEL.filtered.vcf"

BIALLELIC_SNP="${prefix}.biallelic.SNP.raw.vcf"
BIALLELIC_1BP_INDEL="${prefix}.biallelic.1bpINDEL.raw.vcf"
BIALLELIC_gt1BP_INDEL="${prefix}.biallelic.gt1bpINDEL.raw.vcf"

FILTERED_BIALLELIC_SNP="${prefix}.biallelic.SNP.filtered.vcf"
FILTERED_BIALLELIC_1BP_INDEL="${prefix}.biallelic.1bpINDEL.filtered.vcf"
FILTERED_BIALLELIC_gt1BP_INDEL="${prefix}.biallelic.gt1bpINDEL.filtered.vcf"


# select multiallelic sites
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${IN_VCF}" \
    --restrict-alleles-to MULTIALLELIC \
    --output "${MULTIALLELIC}"
# split multiallelic sites
python2 decompose_multiallelic_sites.py \
    "${MULTIALLELIC}" -o "${DECOMPOSED}"

# select multiallelic SNPs
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${DECOMPOSED}" \
    --select-type-to-include SNP \
    --output "${DECOMPOSED_SNP}"
# select multiallelic 1bp indels
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${DECOMPOSED}" \
    --select-type-to-include INDEL \
    --max-indel-size 1 \
    --output "${DECOMPOSED_1BP_INDEL}"
# select multiallelic >1bp indels
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${DECOMPOSED}" \
    --select-type-to-include INDEL \
    --min-indel-size 2 \
    --output "${DECOMPOSED_gt1BP_INDEL}"


# filter multiallelic SNPs
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${DECOMPOSED_SNP}" \
    --output "${FILTERED_DECOMPOSED_SNP}" \
    --filter-name "QDlt2" \
    --filter-expression "AS_QD < 2.0"
# filter multiallelic 1bp indels
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${DECOMPOSED_1BP_INDEL}" \
    --output "${FILTERED_DECOMPOSED_1BP_INDEL}" \
    --filter-name "QDlt5" \
    --filter-expression "AS_QD < 5.0"
# filter multiallelic >1bp indels
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${DECOMPOSED_gt1BP_INDEL}" \
    --output "${FILTERED_DECOMPOSED_gt1BP_INDEL}" \
    --filter-name "QDlt2" \
    --filter-expression "AS_QD < 2.0"

# select biallelic SNPs
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${IN_VCF}" \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    --output "${BIALLELIC_SNP}"
# select biallelic 1bp indels
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${IN_VCF}" \
    --select-type-to-include INDEL \
    --restrict-alleles-to BIALLELIC \
    --max-indel-size 1 \
    --output "${BIALLELIC_1BP_INDEL}"
# select biallelic >1bp indels
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${IN_VCF}" \
    --select-type-to-include INDEL \
    --restrict-alleles-to BIALLELIC \
    --min-indel-size 2 \
    --output "${BIALLELIC_gt1BP_INDEL}"


# filter biallelic SNPs
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${BIALLELIC_SNP}" \
    --output "${FILTERED_BIALLELIC_SNP}" \
    --filter-name "QDlt2" \
    --filter-expression "QD < 2.0"
# filter biallelic 1bp indels
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${BIALLELIC_1BP_INDEL}" \
    --output "${FILTERED_BIALLELIC_1BP_INDEL}" \
    --filter-name "QDlt5" \
    --filter-expression "QD < 5.0"
# filter biallelic >1bp indels
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${BIALLELIC_gt1BP_INDEL}" \
    --output "${FILTERED_BIALLELIC_gt1BP_INDEL}" \
    --filter-name "QDlt2" \
    --filter-expression "QD < 2.0"

# merge files
$GATK MergeVcfs \
    --INPUT "${FILTERED_DECOMPOSED_SNP}" \
    --INPUT "${FILTERED_DECOMPOSED_1BP_INDEL}" \
    --INPUT "${FILTERED_DECOMPOSED_gt1BP_INDEL}" \
    --INPUT "${FILTERED_BIALLELIC_SNP}" \
    --INPUT "${FILTERED_BIALLELIC_1BP_INDEL}" \
    --INPUT "${FILTERED_BIALLELIC_gt1BP_INDEL}" \
    --OUTPUT "${OUT_VCF}"