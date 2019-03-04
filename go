#!/bin/bash

java -Dconfig.file=resources/google.conf -jar bin/cromwell-36.1.jar run wdl/basic_stats/basic_stats.wdl -o resources/gcs_workflow_options.json -i data/smalltest.json 
