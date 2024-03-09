#!/bin/bash

set -euxo pipefail

# synopsis:
# bash /opt/localize_files.sh \
#      object_json \
#      $(pwd) \
#      output_tsv_name

# for localizing files specified in a JSON file,
# formated as follows:
# {
#     "contents":
#     [
#         {
#             "left": "test.pdf",
#             "right": "gs://my_bucket/my_prefix/my_file.txt"
#         }
#     ]
# }

# outputs a 2-col tsv, 
# where 1st col is the desired file name and 
#       2nd col is the local file path
input_json=$1
localize_dir=$2
export localize_dir=${localize_dir%/}
export output_tsv_name=$3

function localize() {
    gcloud storage cp "$1" "${localize_dir}/$2"
    echo -e "$2\t${localize_dir}/$2" >> "${output_tsv_name}"
}


jq --raw-output '.contents[] | .left'  "${input_json}" > tmp.left.txt
jq --raw-output '.contents[] | .right' "${input_json}" > tmp.right.txt
paste tmp.left.txt tmp.right.txt > tmp.2-col.tsv

while IFS=$'\t' read -r -a line
do
    localize "${line[1]}" "${line[0]}"
done < tmp.2-col.tsv

rm tmp.left.txt tmp.right.txt tmp.2-col.tsv
