#!/bin/bash

set -euo pipefail


# reads in metadata, then calls out to create_input_json_for_lrwgs_pipeline.py to produce individual json files according to metadata
export INPUT_PREFIX_EAP="gs://broad-dsde-methods-long-reads-incoming/PBEAP/"
export INPUT_PREFIX_GENENTECH="gs://broad-dsde-methods-long-reads-incoming/Genentech/"
export INPUT_PREFIX_RGP="gs://broad-dsde-methods-long-reads-incoming/RGP/"
export INPUT_PREFIX_VERILY="gs://broad-dsde-methods-long-reads-incoming/Verily/"
export OUTPUT_PREFIX="data/LRWholeGenomeSingleSample/"

parse_and_execute() {
    
    if [[ "eap" == "$2" ]]; then
        INPUT_PREFIX=${INPUT_PREFIX_EAP}
    elif [[ "genentech" == "$2" ]]; then
        INPUT_PREFIX=${INPUT_PREFIX_GENENTECH}
    elif [[ "rgp" == "$2" ]]; then
        INPUT_PREFIX=${INPUT_PREFIX_RGP}
    elif [[ "verily" == "$2" ]]; then
        INPUT_PREFIX=${INPUT_PREFIX_VERILY}
    else
        echo "UNRECOGNIZED MODE: $2."
    fi

    while read -r sp op arr
    do
        if [[ "${sp}" == "SAMPLE_NAME" ]]; then
            continue;
        fi
        
        OUTPUT="${OUTPUT_PREFIX}${op}"
        GCS_PATH=$(echo "${arr}" | sed 's@, @'" ${INPUT_PREFIX}"'@g')
        GCS_PATH="${INPUT_PREFIX}${GCS_PATH}"

        python3 scripts/create_input_json_for_lrwgs_pipeline.py \
          --SM "${sp}" \
          ${GCS_PATH} \
          > "${OUTPUT}"
    done < "${1}"
}
export -f parse_and_execute

parse_and_execute data/LRWholeGenomeSingleSample/eap.tsv eap
parse_and_execute data/LRWholeGenomeSingleSample/genentech.tsv genentech
parse_and_execute data/LRWholeGenomeSingleSample/rgp.tsv rgp
parse_and_execute data/LRWholeGenomeSingleSample/verily.tsv verily
