#!/bin/bash

set -euxo pipefail

export READS=$1
export OVERLAP=$2
export DRAFT=$3

shift 3
export OPTIONS="$@"

racon \
    ${OPTIONS} \
    ${READS} \
    ${OVERLAP} \
    ${DRAFT}
