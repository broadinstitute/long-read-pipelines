#!/usr/bin/python

"""accept_test_results.py

This script copies the results of a CI test ('actual' results) to the
appropriate acceptance test directory ('expected' results).  It is
assumed that the actual directory and the expected directory names
relate as follows:

gs://{actual_bucket}/{workflow_name}/{output...}
gs://{expected_bucket}/test_data/{workflow_name}.test/output_data/{output...}

By default, this script executes in dry-run mode; no files are copied, so
as to provide information on what would be copied.  To actually copy files,
supply the -r flag.

The following is an example command:

python3 scripts/ci/accept_test_results.py \
    -r \
    gs://broad-dsp-lrma-ci/ONT10xSingleFlowcell/ \
    gs://broad-dsp-lrma-ci-resources/test_data

This script calls out to `gsutil rsync`.  Thus, gsutil must be installed
and configured properly.
"""

import argparse
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description='Copy run results to acceptance test dir', prog='accept_test_results')
    parser.add_argument('-r', '--run', action='store_true', help="Run the rsync command")
    parser.add_argument('gcs_act_path', type=str, help="GCS root path to run results")
    parser.add_argument('gcs_exp_path', type=str, help="GCS root path to test results (excluding test name)")
    args = parser.parse_args()

    gcs_act_bucket = re.split("/", re.sub("^gs://", "", args.gcs_act_path))[0]
    gcs_act_prefix = re.sub("^/|/$", "", re.sub(gcs_act_bucket, "", re.sub("^gs://", "", args.gcs_act_path)))

    gcs_exp_dir = f'gs://{gcs_act_bucket}/{gcs_act_prefix}'
    gcs_act_dir = f'gs://{re.sub("^gs://|/$", "", args.gcs_exp_path)}/{gcs_act_prefix}.test/output_data'

    subprocess.run(f'gsutil -m rsync -Cr {"" if args.run else "-n"} {gcs_exp_dir}/ {gcs_act_dir}/', shell=True)


if __name__ == "__main__":
    main()
