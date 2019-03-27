#!/bin/bash

for c in $(gsutil ls "gs://broad-dsde-methods-kiran/workflow_output/CorrectAndAlignWorkflow/" | xargs -L1 basename)
do
    echo $c
    cromshell metadata $c > results/costs/$c.metadata.txt
    python2 scripts/calculate_cost.py -m results/costs/$c.metadata.txt | column -t > results/costs/$c.costs.txt
done
