#!/bin/bash

set -eu

export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

