#!/usr/bin/env python3

import requests
from dmpy import DistributedMake, get_dm_arg_parser
from os.path import basename, dirname
from google.cloud import storage

args = get_dm_arg_parser().parse_args()
m = DistributedMake(args_object=args)
out_dir = 'results/metadata'

storage_client = storage.Client()
bucket = storage_client.get_bucket("broad-dsde-methods-kiran")
blobs = bucket.list_blobs(prefix="workflow_output/CorrectAndAlignWorkflow/", delimiter="/")

for blob in blobs:  # has the side effect of populating blobs.prefixes
    continue

for prefix in blobs.prefixes:
    cromwell_hash = basename(dirname(prefix))

    r_status = requests.get(f'https://cromwell-v36.dsde-methods.broadinstitute.org/api/workflows/v1/{cromwell_hash}/status')
    status = r_status.json()['status']

    if status != 'Running':
        metadata = f'{out_dir}/{cromwell_hash}/{cromwell_hash}.json'
        m.add(metadata, None, f'curl -X GET "https://cromwell-v36.dsde-methods.broadinstitute.org/api/workflows/v1/{cromwell_hash}/metadata?expandSubWorkflows=false" -H "accept: application/json" -s > {metadata}')

        logs_dir = f'{out_dir}/{cromwell_hash}/mlogs/'
        logs_done = f'{logs_dir}/.done'
        m.add(logs_done, metadata, f'python ./scripts/download_monitoring_logs.py {metadata} {logs_dir} && touch {logs_done}')

m.execute()
