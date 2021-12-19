#!/bin/bash

set -eu

for json in *.slim-metadata.json; do
    prefix=$(echo "${json}" | awk -F '.' '{print $1}')
    echo "${prefix}"
    grep -F "id" "${json}" | ggrep -oP '[a-f0-9-]{36}' | sort | uniq > temp.all.txt
    comm -13 "${prefix}.workflow-ids.txt" temp.all.txt > "${prefix}.workflow-ids-misssed.txt"
    rm temp.all.txt
done
