#!/bin/bash

set -eu

executable="/Users/shuang/Projects/work_repos/cromwell-task-monitor-bq/metadata/submit/cromwell_metadata_bq"
if [[ ! -f  "${executable}" ]]; then echo "Cannot find uploader, quite" && exit 1; fi

export CROMWELL_BASEURL="http://34.73.109.5:8000"
export GCP_PROJECT="broad-dsp-lrma"
export DATASET_ID=cromwell_monitoring

echo "===================="
echo "ONT_large"
echo "get slim-metadata"
cromshell slim-metadata "ONT_large" > "ONT_large.slim-metadata.json"
echo "get WID and sub-WID using pure grep"
grep -F "id" "ONT_large.slim-metadata.json" | ggrep -oP '[a-f0-9-]{36}' | sort | uniq > "ONT_large.workflow-ids.txt"
cat "ONT_large.workflow-ids.txt"
echo
echo "  uploading metadata to BQ:${GCP_PROJECT}/${DATASET_ID} ......"
"${executable}" < "ONT_large.workflow-ids.txt"
echo "DONE"
echo "===================="
