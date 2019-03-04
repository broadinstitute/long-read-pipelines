#!/usr/bin/env python3

import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
import json
import subprocess
from utils import parse_workflow_ids

p = argparse.ArgumentParser()
p.add_argument("-p", "--print", action="store_true")
p.add_argument("-o", dest="print", action="store_true")
p.add_argument("workflow_id", nargs="+")
args = p.parse_args()

args.workflow_id = parse_workflow_ids(args.workflow_id)

for workflow_id in args.workflow_id:
    print("="*160)
    print("="*50 + " "*60 + "="*50)
    print("="*50 + f"      WorkflowID: {workflow_id}      " + "="*50)
    print("="*50 + " "*60 + "="*50)
    print("="*160)
    logs = subprocess.check_output(f"cromshell logs {workflow_id}", shell=True)
    logs_json = json.loads(logs)
    for call, shards in sorted(logs_json.get("calls", {}).items()):

        if len(shards) == 0:
            print(f"=== {call} - no shards started")
        elif len(shards) == 1:
            log_path = shards[0]['backendLogs']['log']
            print(f"=== {call}: {log_path}")
        else:
            print(f"=== {call}:")
            for shard in shards:
                attempt_string = f" -- attempt[{shard['attempt']}]" if shard['attempt'] > 1 else ""
                log_path = f"{shard['backendLogs']['log']}"
                print(f"----- shard {shard['shardIndex']}: {log_path}")

        if args.print:
            os.system(f"{'gsutil ' if log_path.startswith('gs') else ''} cat {log_path} | grep -v 'unknown operand'")
            print("------")
            print(f"Log directory:   {os.path.dirname(log_path)}")
            os.system(f"{'gsutil ' if log_path.startswith('gs') else ''} ls {os.path.dirname(log_path)}")
