#!/usr/bin/python

from google.cloud import storage
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


def run_curl_cmd(curl_cmd, shell=False):
    j = {}

    for i in range(3):
        if shell:
            out = subprocess.Popen(curl_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        else:
            out = subprocess.Popen(curl_cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        stdout, stderr = out.communicate()

        if not stdout:
            time.sleep(30)
            exit(1)
        else:
            return json.loads(stdout)

    if not j:
        print_failure(f"No valid response from '{server_url}' for command '{curl_cmd}' after three tries.")
        exit(1)

    return j


def submit_job(wdl, input_json, options, dependencies):
    curl_cmd = f'curl -s -F workflowSource=@{wdl} -F workflowInputs=@{input_json} -F workflowOptions=@{options} -F workflowDependencies=@{dependencies} {server_url}/api/workflows/v1'
    return run_curl_cmd(curl_cmd)


def get_job_status(id):
    curl_cmd = f'curl -s {server_url}/api/workflows/v1/{id}/status'
    return run_curl_cmd(curl_cmd)


def get_job_failure_metadata(id):
    curl_cmd = f'curl -s "{server_url}/api/workflows/v1/{id}/metadata?excludeKey=submittedFiles&expandSubWorkflows=true" | jq ".failures"'

    return run_curl_cmd(curl_cmd, shell=True)


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


def list_blobs(bucket_name, prefix, delimiter=None):
    storage_client = storage.Client()
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix, delimiter=delimiter)

    return blobs


def upload_metadata(test, id):
    md_cmd = f'curl -s "{server_url}/api/workflows/v1/{id}/metadata?excludeKey=submittedFiles&expandSubWorkflows=true" | jq . > md.txt'
    subprocess.run(md_cmd, shell=True)

    storage_client = storage.Client()
    bucket = storage_client.bucket('broad-dsp-lrma-ci')
    blob = bucket.blob(f'metadata/{test}/{id}.metadata.txt')

    blob.upload_from_filename('md.txt')

    return f'gs://broad-dsp-lrma-ci/metadata/{test}/{id}.metadata.txt'


def find_outputs(input_json, exp_bucket='broad-dsp-lrma-ci-resources', act_bucket='broad-dsp-lrma-ci'):
    b = os.path.basename(input_json).replace(".json", "")
    outs = {}

    expected_blobs = list_blobs(exp_bucket, f'test_data/{b}/output_data')
    for blob in expected_blobs:
        bn = os.path.basename(blob.name)

        outs[bn] = {'exp': blob.md5_hash, 'exp_path': f'gs://{exp_bucket}/{blob.name}', 'act': None, 'act_path': None}

    with open(input_json) as jf:
        for l in jf:
            if f'gs://{act_bucket}/' in l:
                m = re.split(":\\s+", re.sub("[\",]+", "", l.strip()))
                if len(m) > 1:
                    n = re.sub("/$", "", m[-1].replace(f'gs://{act_bucket}/', ""))

                    blobs = list_blobs(act_bucket, n)
                    for blob in blobs:
                        bn = os.path.basename(blob.name)

                        if bn not in outs:
                            print_warning(f"Found an actual output file {bn} that is not in the expected directory")
                            continue

                        outs[bn]['act'] = blob.md5_hash
                        outs[bn]['act_path'] = f'gs://{act_bucket}/{blob.name}'

    return outs


def compare_contents(exp_path, act_path):
    storage_client = storage.Client()

    fn, ext = os.path.splitext(exp_path)

    exp = f'exp{ext}'
    act = f'act{ext}'

    with open(exp, "wb") as exp_obj:
        storage_client.download_blob_to_file(exp_path, exp_obj)

    with open(act, "wb") as act_obj:
        storage_client.download_blob_to_file(act_path, act_obj)

    if ext == '.fastq' or exp_path.endswith('.fastq.gz') or exp_path.endswith('.fq.gz') or ext == '.fasta' or exp_path.endswith('.fasta.gz') or exp_path.endswith('.fa.gz'):
        r = subprocess.run(f'mash dist -t {exp} {act} 2>/dev/null | grep -v "query" | awk "{{ exit $2 }}"', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
        if ext == '.bam':
            subprocess.run(f'samtools view {exp} | sort > exp.tmp', shell=True)
            subprocess.run(f'samtools view {act} | sort > act.tmp', shell=True)
        elif ext == '.gz':
            subprocess.run(f'zcat {exp} | grep -v -e fileDate > exp.tmp', shell=True)
            subprocess.run(f'zcat {act} | grep -v -e fileDate > act.tmp', shell=True)
        elif ext == '.pdf':
            subprocess.run(f'pdftotext {exp} > exp.tmp', shell=True)
            subprocess.run(f'pdftotext {act} > act.tmp', shell=True)
        else:
            print_warning(f'Unknown file extension {ext} for file {exp_path} and {act_path}')
            return 1

        r = subprocess.run(f'diff exp.tmp act.tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        os.remove('exp.tmp')
        os.remove('act.tmp')

    os.remove(exp)
    os.remove(act)

    if r is not None:
        if r.stdout != b'' or r.stderr != b'':
            print_warning(f'comparing "{exp_path}" vs "{act_path}"')
            print_warning(r.stdout.decode('utf-8'))
            print_warning(r.stderr.decode('utf-8'))
            print_warning(r.returncode)

        return r.returncode

    return 1


def compare_outputs(test, outs):
    print_info(f'{test}')

    num_mismatch = 0
    for b in outs:
        if not b.endswith(".png") and not b.endswith("sequencing_summary.txt") and outs[b]['exp'] != outs[b]['act']:
            if outs[b]['act_path'] is None or compare_contents(outs[b]['exp_path'], outs[b]['act_path']) != 0:
                print_info(f'- {b} versions are different:')
                print_failure(f"    exp: ({outs[b]['exp']}) {outs[b]['exp_path']}")
                print_failure(f"    act: ({outs[b]['act']}) {outs[b]['act_path']}")
                num_mismatch += 1

    return num_mismatch


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
    old_num_finished = -1

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

            if old_num_finished != num_finished:
                print_info(f'Running {len(jobs)} tests, {num_finished} tests complete. {num_succeeded} succeeded, {num_failed} failed.')

            old_num_finished = num_finished
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
            mdpath = upload_metadata(test, jobs[test]["id"])

            print_success(f"{test}: Workflow {jobs[test]['status']} ({diff.total_seconds()}s -- {str(diff)})")
            print_success(f"{test}: Workflow metadata uploaded to {mdpath}")

            outs = find_outputs(input[test])
            num_mismatch = compare_outputs(test, outs)

            if num_mismatch == 0:
                print_success(f"{test}: {len(outs)} files checked, {num_mismatch} failures")
            else:
                print_failure(f"{test}: {len(outs)} files checked, {num_mismatch} failures")

                ret = 1
        else:
            mdpath = upload_metadata(test, jobs[test]["id"])

            print_failure(f"{test}: Workflow {jobs[test]['status']} ({diff.total_seconds()}s -- {str(diff)})")
            print_failure(f"{test}: Workflow metadata uploaded to {mdpath}")
            print_failure(f"{test}: Workflow failure messages:\n{json.dumps(get_job_failure_metadata(jobs[test]['id']), sort_keys=True, indent=4)}")

            ret = 1


if ret == 0:
    print_success('ALL TESTS SUCCEEDED.')
else:
    print_failure('SOME TESTS FAILED.')

exit(ret)
