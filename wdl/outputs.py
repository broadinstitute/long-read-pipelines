#!/usr/bin/env python3

import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
import json
import subprocess
from utils import parse_workflow_ids

p = argparse.ArgumentParser()
p.add_argument("-d", "--download", help="Download output to this path")
p.add_argument("workflow_id", nargs="+")
args = p.parse_args()

args.workflow_id = parse_workflow_ids(args.workflow_id)

if args.download and not os.path.isdir(args.download):
    p.exit(f"Directory doesn't exist: {args.download}")

for workflow_id in args.workflow_id:
    print("="*160)
    print("="*50 + " "*60 + "="*50)
    print("="*50 + f"      WorkflowID: {workflow_id}      " + "="*50)
    print("="*50 + " "*60 + "="*50)
    print("="*160)
    metadata = subprocess.check_output(f"cromshell metadata {workflow_id}", shell=True)
    metadata_json = json.loads(metadata)
    for call, outputs in metadata_json.get("outputs", {}).items():
        if isinstance(outputs, str):
            if args.download:
                if not any([outputs.endswith(suffix) for suffix in [".bam", ".fa", ".fasta"]]):
                    os.system(f"gsutil -m cp {outputs} {args.download}")
                else:
                    print(f"=== {call}: Skipping download of {outputs}")

        else:
            print(f"=== {call}:")
            for output in outputs[0:3]:
                print(f"----- output {output}")
            if len(outputs) > 3:
                print("...")