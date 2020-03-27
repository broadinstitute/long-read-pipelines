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


def print_test_failure(test, message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'\033[1;31;40m[{ts} FAIL] {test}: {message}\033[0m')


def print_test_success(test, message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'\033[1;32;40m[{ts} PASS] {test}: {message}\033[0m')


def print_test_warning(test, message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'\033[1;33;40m[{ts} WARN] {test}: {message}\033[0m')


def print_test_progress(test, message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'[{ts} INFO] {test}: {message}')


def print_info(message):
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f'[{ts} INFO] {message}')


def list_test_inputs():
    return glob.glob("test/test_json/*.json")


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
        j = json.loads(stdout)

        if 'id' in j:
            return j

        time.sleep(30)

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
        if jobs[test]["status"] == "Running":
            running = True

    return running


print_info(f'Cromwell server: {server_url}')

input_jsons = list_test_inputs()
print_info(f'Found {len(input_jsons)} tests.')

ret = 0

jobs = {}
for input_json in input_jsons:
    test = os.path.basename(input_json)
    p = test.split(".")

    wdl_name = f'{p[0]}.wdl'
    wdl_path = find_wdl_path(wdl_name)

    if wdl_path is None:
        print_test_warning(test, f'Requested WDL does not exist.')
    else:
        if 'TestCromwell' in wdl_path:
            j = submit_job(wdl_path, input_json, 'resources/workflow_options/ci.json', 'wdl/lr_wdls.zip')

            print_test_progress(test, f'{j["id"]}, {j["status"]}')
            jobs[test] = j

if len(jobs) > 0:
    num_finished = len(jobs)
    num_failed = 0
    num_succeeded = 0

    while True:
        time.sleep(60)
        jobs = update_status(jobs)

        if jobs_are_running(jobs):
            num_finished = len(jobs)
            num_failed = 0
            num_succeeded = 0

            for test in jobs:
                if jobs[test]['status'] == 'Running':
                    num_finished = num_finished - 1
                elif jobs[test]['status'] == 'Failed':
                    num_failed = num_failed + 1
                elif jobs[test]['status'] == 'Succeeded':
                    num_succeeded = num_succeeded + 1

            print_info(f'Running {len(jobs)} tests. {num_finished} tests complete, {num_succeeded} succeeded, {num_failed} failed.')
        else:
            break

    print_info(f'Ran {len(jobs)} tests. {num_finished} tests complete, {num_succeeded} succeeded, {num_failed} failed.')

    for test in jobs:
        if jobs[test]['status'] == 'Succeeded':
            print_test_success(test, jobs[test]['status'])
        else:
            print_test_failure(test, jobs[test]['status'])
            ret = 1

if ret == 0:
    print_info('ALL TESTS SUCCEEDED.')

exit(ret)
