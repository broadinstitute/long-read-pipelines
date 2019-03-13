#!/bin/bash

set -euxo pipefail

java -jar bin/womtool-36.1.jar validate wdl/correct_and_align/correct_and_align.wdl

#cromshell submit wdl/correct_and_align/correct_and_align.wdl data/smalltest.json resources/gcs_workflow_options.json
cromshell submit wdl/correct_and_align/correct_and_align.wdl data/ecoli.json resources/gcs_workflow_options.json
