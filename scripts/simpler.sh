#!/bin/bash

set -euo pipefail


# reads in metadata, then calls out to create_input_json_for_lrwgs_pipeline.py to produce individual json files according to metadata
export INPUT_PREFIX_EAP=" gs://broad-dsde-methods-kiran/pb_eap/"
export INPUT_PREFIX_INC=" gs://broad-dsde-methods-long-reads/incoming/"
export OUTPUT_PREFIX="../data/LRWholeGenomeSingleSample/"

parse_and_execute() {
    
    if [[ "eap" == "$2" ]]; then
        INPUT_PREFX=${INPUT_PREFIX_EAP}
    elif [[ "incoming" == "$2" ]]; then
        INPUT_PREFX=${INPUT_PREFIX_INC}
    else
        echo "UNRECOGANIZED MODE: $2."
    fi

    while read -r sp op arr
    do
        if [[ "${sp}" == "SAMPLE_NAME" ]]; then
            continue;
        fi
        
        OUTPUT="${OUTPUT_PREFIX}${op}"
        GCS_PATH=$(echo "${arr}" | sed 's@, @'" ${INPUT_PREFX}"'@g')
        GCS_PATH="${INPUT_PREFX}${GCS_PATH}"

        python3 create_input_json_for_lrwgs_pipeline.py \
          --SM "${sp}" \
          "${GCS_PATH}" \
          > "${OUTPUT}"
    done < "${1}"
}
export -f parse_and_execute

parse_and_execute ../data/LRWholeGenomeSingleSample/eap.tsv eap
parse_and_execute ../data/LRWholeGenomeSingleSample/incoming.tsv incoming
