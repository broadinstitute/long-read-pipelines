#!/usr/bin/env python3

import os

os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
import datetime
import json
import re
from pprint import pformat
from subprocess import check_output, DEVNULL
from utils import parse_workflow_ids

def parse_date(isoformat_date_string):
    if isoformat_date_string == "now":
        return datetime.datetime.now()

    dt, _, us = isoformat_date_string.rstrip('Z').partition(".")
    #us = int(us.rstrip("Z"))
    return datetime.datetime.strptime(dt, "%Y-%m-%dT%H:%M:%S") # + datetime.timedelta(microseconds=us)

def get_text_file_length(path, filters=()):
    if path.startswith("gs://"):
        command = f"gsutil cat {path}"
    elif os.path.isfile(path):
        command = f"cat {path}"
    else:
        print("ERROR: invalid path: {path}")
        return 0

    if path.endswith("gz"):
        command += " | gunzip -c -"

    command += f" | {' | '.join(['grep -v ^#'] + list(filters))}"
    command += " | wc -l"

    #print(command)
    result = check_output(command, shell=True)

    return int(result)


def print_paths_of_interest(key, path, fast=False):

    if isinstance(path, str) and any(suffix in path for suffix in [".bed", ".vcf", ".bam", ".bai", ".log"] if not any(suffix in path for suffix in [".tbi"])):
        summary_str = f"{key}: {path}"
        if any(suffix in path for suffix in [".bed", ".vcf"]):
            num_calls = ""
            if not fast:
                if ".vcf" in path:
                    num_calls = f", calls: {get_text_file_length(path, filters=['grep -v 0/0:']):,}"  # 'bcftools filter -e GT=RR'
                summary_str += f" ({get_text_file_length(path):,} records{num_calls})"
        print(summary_str)

def summarize_call(call_name, call_json):
    # keys:
    #      'shardIndex', 'attempt', 'jobId', 'backendLabels', 'labels',
    #      'start', 'end',
    #      'inputs', 'outputs',
    #      'stdout', 'stderr', 'returnCode',
    #      'callRoot', 'executionEvents',
    #      'executionStatus', 'backendStatus',  'backendLogs',
    #      'jes', 'runtimeAttributes', 'callCaching', 'dockerImageUsed', 'backend',
    overall_timedeltas = []
    main_command_timedeltas = []
    localization_timedeltas = []
    log_files = []
    for i, scatter_json in enumerate(call_json):
        if scatter_json.get('backendLogs'):
            log_files.append(f"{scatter_json.get('backendLogs', {}).get('log', '')}")
        if scatter_json.get('backendStatus') == "Failed" or scatter_json.get('executionStatus') == "Failed":
            print(f"                      ==> ERROR: {call_name} shard #{scatter_json['shardIndex']} failed:\n\n"
                f"{pformat(scatter_json.get('failures', ''))}\n\n")
            continue
        if 'start' in scatter_json:
            overall_timedeltas.append(parse_date(scatter_json.get('end', 'now')) - parse_date(scatter_json['start']))
        for execution_event in scatter_json.get('executionEvents', []):
            duration_timedetla = parse_date(execution_event['endTime']) - parse_date(execution_event['startTime'])
            if duration_timedetla.total_seconds() < 5*60:
                continue

            is_main_command = re.search('Started running "/bin/bash /cromwell_root/script"', execution_event['description'])
            if is_main_command:
                main_command_timedeltas.append(duration_timedetla)

            is_localization = re.search('Started.+gsutil[ ]+cp[ ]+(gs://[^ ]+)', execution_event['description'])
            if is_localization:
                remote_file_path = is_localization.group(1)
                localization_timedeltas.append(duration_timedetla)

            #if is_localization or is_main_command:
            #    print(f"{parse_date(execution_event['endTime']) - parse_date(execution_event['startTime'])}     {execution_event['description'][:170]} is_main_command={is_main_command}, is_localization={is_localization}")

    time_string = ""
    for timedeltas, label in [(overall_timedeltas, 'run'), (localization_timedeltas, 'localization'), (main_command_timedeltas, 'command')]:
        if not timedeltas:
            continue
        timedeltas_sum = sum(timedeltas, datetime.timedelta(0))
        time_string += f" -  mean {label} time: {str(timedeltas_sum / len(timedeltas)).split('.')[0]}  "
        if len(call_json) > 1:
            time_string += f"[{min(timedeltas)} : {max(timedeltas)}],"
            time_string += f" sum: {timedeltas_sum.total_seconds()/3600.:0.1f}h "

    print(f"{len(call_json):2d} x {call_name:39s}{time_string}")
    print(("\n".join(map(lambda l: 25*" "+l, log_files[0:3]))) if log_files else "")


def summarize_run(workflow_id, fast=False):
    command = f"cromshell metadata {workflow_id}"
    print(f"{command}\n")
    result = check_output(command, stderr=DEVNULL, shell=True)

    try:
        data = json.loads(result)
    except Exception as e:
        os.system(f"ERROR: unable to parse: {e}\n{result}")

    for call_key, call_json in data['calls'].items():
        summarize_call(call_key, call_json)

    print("")
    print(59*" " + f"- total wall clock time: {parse_date(data.get('end', 'now')) - parse_date(data['start'])}")
    print(59*" " + f"- status: {data['status']}")

    print("")
    for key, value in sorted(list(data["inputs"].items()) + list(data["outputs"].items()), key=lambda x: str(x[1]).partition('.')[-1]):
        print_paths_of_interest(key, value, fast=fast)

    print("")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-f", "--fast", action="store_true")
    p.add_argument("workflow_id", nargs="+")
    args = p.parse_args()

    args.workflow_id = parse_workflow_ids(args.workflow_id)

    for wid in args.workflow_id:
        print("-------------------")
        summarize_run(wid, fast=args.fast)
