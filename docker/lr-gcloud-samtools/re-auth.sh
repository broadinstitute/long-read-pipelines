#!/bin/bash

set -eu

export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
export GCS_REQUESTER_PAYS_PROJECT=$(gcloud config get-value project)
