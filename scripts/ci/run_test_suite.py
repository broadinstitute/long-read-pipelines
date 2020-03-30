#!/usr/bin/python

import glob
import os
import re
import subprocess
import json
import time
import datetime
import sys


server_url = "http://localhost:8000"
if len(sys.argv) == 2:
    server_url = sys.argv[1]


def print_failure(message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'\033[1;31;40m[{ts} FAIL] {message}\033[0m')


def print_success(message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'\033[1;32;40m[{ts} PASS] {message}\033[0m')


def print_warning(message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'\033[1;33;40m[{ts} WARN] {message}\033[0m')


def print_info(message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'[{ts} INFO] {message}')


def list_test_inputs():
    return glob.glob("test/test_json/*.json")


def list_disabled_tests():
    disabled_tests = set()

    disabled_tests_path = "test/test_json/all_disabled_tests.txt"
    if os.path.exists(disabled_tests_path):
        with open(disabled_tests_path) as f:
            for line in f:
                disabled_tests.add(line.strip())

    return disabled_tests


def find_wdl_path(wdl):
    for (root, dirs, files) in os.walk("wdl/"):
        for file in files:
            my_wdl = re.sub("/+", "/", f'{root}/{file}')
            if my_wdl.endswith(wdl):
                return my_wdl

    return None


def run_curl_cmd(curl_cmd):
    j = {}

    for i in range(3):
        out = subprocess.Popen(curl_cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = out.communicate()

        if not stdout:
            print_failure(f"Unable to dispatch command '{curl_cmd}' to '{server_url}'")
            exit(1)

        j = json.loads(stdout)

        if 'id' in j:
            return j

        time.sleep(30)

    if 'id' not in j:
        print_failure(f"No valid response from '{server_url}' for command '{curl_cmd}' after three tries.")
        exit(1)

    return j


def submit_job(wdl, input_json, options, dependencies):
    curl_cmd = f'curl -s -F workflowSource=@{wdl} -F workflowInputs=@{input_json} -F workflowOptions=@{options} -F workflowDependencies=@{dependencies} {server_url}/api/workflows/v1'
    return run_curl_cmd(curl_cmd)


def get_job_status(id):
    curl_cmd = f'curl -s {server_url}/api/workflows/v1/{id}/status'
    return run_curl_cmd(curl_cmd)


def update_status(jobs):
    new_jobs = {}
    for test in jobs:
        new_jobs[test] = get_job_status(jobs[test]["id"])

    return new_jobs


def jobs_are_running(jobs):
    running = False

    for test in jobs:
        if jobs[test]["status"] == "Running" or jobs[test]["status"] == "Submitted":
            running = True

    return running


def find_outputs(input_json):
    b = os.path.basename(input_json).replace(".json", "")
    reference_output_dir = f'gs://broad-dsp-lrma-ci-resources/test_data/{b}/output_data'
    p1 = subprocess.Popen(f'gsutil hash {reference_output_dir}/**'.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p2 = subprocess.Popen(f'paste - - -'.split(' '), stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = p2.communicate()

    outs = {}
    for f in stdout.decode('utf-8').split('\n'):
        g = re.split("\s+", f.strip())

        if len(g) > 1:
            bn = os.path.basename(g[3].replace(":", ""))
            ha = g[len(g) - 1]
            outs[bn] = {'exp': ha, 'exp_path': g[3], 'act': None, 'act_path': None}

    with open(input_json) as jf:
        for l in jf:
            if 'gs://broad-dsp-lrma-ci' in l:
                m = re.split(":\\s+", re.sub("[\",]+", "", l.strip()))
                if len(m) > 1:
                    n = re.sub("/+$", "", m[1])

                    p1 = subprocess.Popen(f'gsutil hash {n}/**'.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    p2 = subprocess.Popen(f'paste - - -'.split(' '), stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    stdout, stderr = p2.communicate()

                    for h in stdout.decode('utf-8').split('\n'):
                        i = re.split("\s+", h.strip())

                        if len(i) > 1:
                            bm = os.path.basename(i[3].replace(":", ""))
                            hb = i[len(i) - 1]

                            if bm in outs:
                                outs[bm]['act'] = hb
                                outs[bm]['act_path'] = i[3]
                            else:
                                outs[bm] = {'exp': None, 'exp_path': None, 'act': hb, 'act_path': i[3]}

    outs.pop('matched', None)

    return outs


def compare_outputs(test, outs):
    ret = 0
    for b in outs:
        if outs[b]['exp'] != outs[b]['act']:
            print_failure(f"{test}: {outs[b]['exp_path']} ({outs[b]['exp']} != {outs[b]['act_path']} ({outs[b]['act']}")
            ret = 1

    return ret


print_info(f'Cromwell server: {server_url}')

# List tests
input_jsons = list_test_inputs()
disabled_tests = list_disabled_tests()
print_info(f'Found {len(input_jsons)} tests, {len(disabled_tests)} disabled.')
for input_json in input_jsons:
    if input_json in disabled_tests:
        print_warning(f'[ ] {input_json}')
    else:
        print_info(f'[*] {input_json}')

# Dispatch tests
jobs = {}
times = {}
input = {}
for input_json in input_jsons:
    if input_json not in disabled_tests:
        test = os.path.basename(input_json)
        p = test.split(".")

        wdl_name = f'{p[0]}.wdl'
        wdl_path = find_wdl_path(wdl_name)

        if wdl_path is None:
            print_warning(f'{test}: Requested WDL does not exist.')
        else:
            j = submit_job(wdl_path, input_json, 'resources/workflow_options/ci.json', 'wdl/lr_wdls.zip')

            print_info(f'{test}: {j["id"]}, {j["status"]}')
            jobs[test] = j
            times[test] = {'start': datetime.datetime.now(), 'stop': None}
            input[test] = input_json

# Monitor tests
ret = 0
if len(jobs) > 0:
    while True:
        time.sleep(60)
        jobs = update_status(jobs)

        if jobs_are_running(jobs):
            num_finished = len(jobs)
            num_failed = 0
            num_succeeded = 0

            for test in jobs:
                if jobs[test]['status'] == 'Running' or jobs[test]['status'] == 'Submitted':
                    num_finished = num_finished - 1
                elif jobs[test]['status'] == 'Failed':
                    num_failed = num_failed + 1
                    if times[test]['stop'] is None:
                        times[test]['stop'] = datetime.datetime.now()
                elif jobs[test]['status'] == 'Succeeded':
                    num_succeeded = num_succeeded + 1
                    if times[test]['stop'] is None:
                        times[test]['stop'] = datetime.datetime.now()

            print_info(f'Running {len(jobs)} tests, {num_finished} tests complete. {num_succeeded} succeeded, {num_failed} failed.')
        else:
            break

    num_finished = len(jobs)
    num_failed = 0
    num_succeeded = 0

    for test in jobs:
        if jobs[test]['status'] == 'Succeeded':
            num_succeeded = num_succeeded + 1
        else:
            num_failed = num_failed + 1

        if times[test]['stop'] is None:
            times[test]['stop'] = datetime.datetime.now()

    print_info(f'Finished {num_finished} tests. {num_succeeded} succeeded, {num_failed} failed.')

    for test in jobs:
        diff = times[test]['stop'] - times[test]['start']
        if jobs[test]['status'] == 'Succeeded':
            print_success(f"{test}: {jobs[test]['status']} ({diff.total_seconds()}s -- {str(diff)})")

            outs = find_outputs(input[test])
            ret = compare_outputs(test, outs)
        else:
            print_failure(f"{test}: {jobs[test]['status']} ({diff.total_seconds()}s -- {str(diff)})")
            ret = 1


if ret == 0:
    print_success('ALL TESTS SUCCEEDED.')
else:
    print_failure('SOME TESTS FAILED.')

exit(ret)
