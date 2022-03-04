#!/bin/bash

set -eu

################################################################################
# A script to collect which dockers are in use and which latest dockers available
################################################################################

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd "${dir}"


echo "COLLECTING DOCKERS IN USE..."
cd ../wdl
rm -f dockers.in_use.tsv
for wdl in $(find . -name "*.wdl"| sed "s|^\./||") ; do
    if ! grep -qE 'docker:\s+\"' "${wdl}"; then continue; fi;
    grep -nE 'docker:\s+\"' "${wdl}" > tmp.0.txt
    awk -F ':' '{print $3"\t"$4"\t"$1}' tmp.0.txt | sed -e 's/^[[:space:]]*//' | sed "s/\"//g" | awk -F '/' '{print $NF}' | sort > tmp.1.txt
    sed -e "s%$%\t$wdl%" tmp.1.txt >> dockers.in_use.tsv
    rm tmp.*.txt
done
echo -e "name\ttag\tline\twdl\n$(sort dockers.in_use.tsv)" > dockers.in_use.sorted.tsv
rm dockers.in_use.tsv

echo "COLLECTING LATEST DOCKERS AVAILABLE..."
cd ../docker
rm -f dockers.latest.tsv
for makefile in $(find . -mindepth 2 -name "Makefile" | sed "s|^\./||") ; do
    name=$(grep -m 1 -F 'TAG1' "${makefile}" | awk -F '/' '{print $NF}' | awk -F ':' '{print $1}')
    tag=$(head -n 1 "${makefile}" | awk -F '=' '{print $NF}' | sed 's% %%g' | awk -F '#' '{print $1}')
    echo -e "${name}\t${tag}" >> dockers.latest.tsv
done
sort dockers.latest.tsv > dockers.latest.sorted.tsv
rm dockers.latest.tsv

echo "DONE. PLEASE CHECKOUT TWO TSV FILES: [dockers.in_use.sorted.tsv, dockers.latest.sorted.tsv]"
cd "${dir}"
mv ../wdl/dockers.in_use.sorted.tsv .
mv ../docker/dockers.latest.sorted.tsv .
