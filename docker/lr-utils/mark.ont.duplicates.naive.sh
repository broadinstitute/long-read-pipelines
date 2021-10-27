#!/bin/bash

set -eu

if [ "$#" -lt 1 ]; then
  echo "I need a BAM file to work on" && exit 1
fi

BAM=$1

filename=$(basename -- "${BAM}")
# extension="${filename##*.}"
prefix="${filename%.*}"

echo "==========================================================="
echo "collecting duplicate information" && date
time \
    samtools view -@ 1 "${BAM}" | \
    awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' | \
    sort | uniq -d \
    > "${BAM}".duplicates.txt
echo "==========================================================="
echo "de-duplicating" && date
time python3 remove_duplicate_ont_aln.py \
    --prefix "${prefix}" \
    --annotations "${BAM}".duplicates.txt \
    "${BAM}"
echo "==========================================================="
echo "DONE"
