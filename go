#!/bin/bash

#java -Dconfig.file=resources/google.conf -jar bin/cromwell-36.1.jar run wdl/basic_stats/basic_stats.wdl -o resources/gcs_workflow_options.json -i data/smalltest.json 

#cromshell submit wdl/basic_stats/basic_stats.wdl data/smalltest.json resources/gcs_workflow_options.json
cromshell submit wdl/correct_and_align/correct_and_align.wdl data/smalltest.json resources/gcs_workflow_options.json
