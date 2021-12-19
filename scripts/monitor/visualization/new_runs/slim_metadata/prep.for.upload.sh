#!/bin/bash

set -eu

submissions=( ASM_small_final ASM_medium_final ASM_large_final CCS_large_96G_16cores CCS_small ONT_small )

for sub in "${submissions[@]}"; do
    echo "===================="
    echo "${sub}"
    echo "get slim-metadata"
    cromshell slim-metadata "${sub}" > "${sub}.slim-metadata.json"
    echo "get WID and sub-WID"
    /Users/shuang/Projects/work_repos/pipelines/scripts/monitor/utilities/get_workflow_subworkflow_ids.py \
        "${sub}.slim-metadata.json" \
        > "${sub}.WID-and-subWID.json"
    echo "parse WIDs into txt"
    ggrep -oP '[a-f0-9-]{36}' "${sub}.WID-and-subWID.json" | sort | uniq > "${sub}.workflow-ids.txt"
    echo "DONE"
    echo "===================="
done

