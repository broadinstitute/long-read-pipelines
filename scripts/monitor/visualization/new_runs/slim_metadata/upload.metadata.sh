#!/bin/bash

set -eu

executable="/Users/shuang/Projects/work_repos/cromwell-task-monitor-bq/metadata/submit/cromwell_metadata_bq"
if [[ ! -f  "${executable}" ]]; then echo "Cannot find uploader, quite" && exit 1; fi

export CROMWELL_BASEURL="http://34.73.109.5:8000"
export GCP_PROJECT="broad-dsp-lrma"
export DATASET_ID=cromwell_monitoring

for ff in *.workflow-ids.txt; do
    echo "===================="
    prefix=$(echo ${ff} | awk -F '.' '{print $1}')
    echo ${prefix}
    echo "  uploading metadata to BQ:${GCP_PROJECT}/${DATASET_ID} ......"
    "${executable}" < "${ff}"
    echo "DONE"
    echo "===================="
done