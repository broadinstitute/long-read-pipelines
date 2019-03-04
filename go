#!/bin/bash

java -Dconfig.file=resources/google.conf -jar bin/cromwell-36.1.jar run wdl/basic_stats/basic_stats.wdl -i data/ecoli.json -o resources/gcs_workflow_options.json
