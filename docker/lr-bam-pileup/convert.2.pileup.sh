#!/bin/bash

set -eu

input_bam=$1
ref_fasta=$2
baq_option=$3
my_bed=$4

prefix=$(echo "${my_bed}" | awk -F '.' '{print $1}')

samtools view -h \
    --region-file "${my_bed}" \
    --write-index \
    -o "${prefix}.bam##idx##${prefix}.bam.bai" \
    "${input_bam}"

samtools mpileup \
    "${baq_option}" \
    -s \
    -q 1 \
    -f "${ref_fasta}" \
    -o "${prefix}.mpileup" \
    "${prefix}.bam" \
    2> "${prefix}.mpileup.err"

rm "${prefix}.bam" "${prefix}.bam.bai"
