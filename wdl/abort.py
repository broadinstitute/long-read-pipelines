#!/usr/bin/env python3

import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
from utils import parse_workflow_ids

p = argparse.ArgumentParser()
p.add_argument("workflow_id", nargs="+")
args = p.parse_args()

args.workflow_id = parse_workflow_ids(args.workflow_id)

for workflow_id in args.workflow_id:
    os.system(f"cromshell abort {workflow_id}")
