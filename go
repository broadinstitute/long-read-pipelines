#!/bin/bash

set -euxo pipefail

java -jar bin/womtool-36.1.jar validate wdl/correct_and_align/correct_and_align.wdl
java -jar bin/womtool-36.1.jar validate wdl/compute_metrics/compute_metrics.wdl

#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/smalltest.json            resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/ecoli.json                resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/CEU_NA12878rep1.json      resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/CEU_NA12878rep2.json      resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/CEU_NA12891.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/CEU_NA12892.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/CHS_HG00512.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/CHS_HG00513.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/CHS_HG00514.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/GWD_HG02982_CLR.json      resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/GWD_HG02983.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/GWD_HG02984.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/YRI_NA19238.json          resources/gcs_workflow_options.json
#sleep 5
#cromshell  submit  wdl/correct_and_align/correct_and_align.wdl  data/YRI_NA19239.json          resources/gcs_workflow_options.json
#sleep 5

cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/smalltest.json            resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/ecoli.json                resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/CEU_NA12878rep1.json      resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/CEU_NA12878rep2.json      resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/CEU_NA12891.json          resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/CEU_NA12892.json          resources/gcs_workflow_options_metrics.json
sleep 5
#cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/CHS_HG00512.json          resources/gcs_workflow_options_metrics.json
#sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/CHS_HG00513.json          resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/CHS_HG00514.json          resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/GWD_HG02982_CLR.json      resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/GWD_HG02983.json          resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/GWD_HG02984.json          resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/YRI_NA19238.json          resources/gcs_workflow_options_metrics.json
sleep 5
cromshell  submit  wdl/compute_metrics/compute_metrics.wdl  data/YRI_NA19239.json          resources/gcs_workflow_options_metrics.json
sleep 5
