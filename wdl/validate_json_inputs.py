#!/usr/bin/env python3

import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
import json
import sys
from utils import parse_workflow_ids

def gsutil_ls(gs_path):
    result = os.system(f"gsutil ls {gs_path}")
    if result != 0:
        print(f"ERROR: Invalid path {gs_path}")
        return 1
    return 0

def validate_input_json(cromwell_input_json_file_path):
    if not os.path.isfile(cromwell_input_json_file_path):
        raise ValueError(f"Invalid input file path: {cromwell_input_json_file_path}")

    error_counter = 0
    with open(cromwell_input_json_file_path) as f:
        for key, value in json.load(f).items():
            if not isinstance(value, str):
                continue
            print(f"==> key: {key}: {value}")
            if value.startswith("gs://"):
                error_counter += gsutil_ls(value)
            elif os.path.isabs(value) and not os.path.exists(value):
                print(f"local file not found: {value}")
                error_counter += 1

    if error_counter:
        raise ValueError(f"{error_counter} errors found")
    else:
        print(f"success: 0 errors found")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workflow_id", nargs="+")
    args = p.parse_args()

    args.workflow_id = parse_workflow_ids(args.workflow_id)

    for cromwell_input_json_file_path in args.workflow_id:
        print("="*100)
        print(f"    validating {cromwell_input_json_file_path}")
        print("="*100)

        try:
            validate_input_json(cromwell_input_json_file_path)
        except ValueError as e:
            sys.exit(f"ERROR: {str(e)}")