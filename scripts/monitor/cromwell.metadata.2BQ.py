#!/usr/bin/env python3

import os
import sys
import errno
import subprocess
import pandas as pd
from shutil import copyfile

HOME = os.getenv('HOME')

###############################
# user-configurable variables #
###############################

gcp_project = "broad-dsde-methods"
work_dir = f"{HOME}/Desktop/cromwell_monitoring"

#####################
# derived variables #
#####################

# throw an error if cromwell-task-monitor-bq is unavailable
try:
    with open("/dev/null", "w") as hide: # pipe output to /dev/null for silence
        subprocess.Popen("cromwell_metadata_bq", stdout=hide, stderr=hide)
except OSError:
    print("cromwell_metadata_bq not found")

# check requested files exist
if ( not os.path.exists(work_dir) ):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            work_dir)

cromshell_database = f"{HOME}/.cromshell/all.workflow.database.tsv"
if ( not os.path.exists(cromshell_database) ):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            cromshell_database)

# Get cromwell server info from cromshell files
with open(HOME + "/.cromshell/cromwell_server.config", "r") as f:
    cromwell_base_url = f.readline().rstrip()

###################################
# Create or Update WorkflowID.txt #
###################################

workflow_id_file = work_dir + '/workflowIds.txt'
workflow_id_tmp_file = work_dir + '/workflowIds_temp.txt'

# Overwrite workflow_id_file if it exists, othewise create a new one
print("\n## Obtaining Wokflow IDs ##")
if (os.path.exists(workflow_id_file) and
    0 < len(open(workflow_id_file).readlines())):

    print("Grab last wokflow ID whose metadata was already pushed")
    last_workflow_id = ""
    with open(workflow_id_file, "r") as f:
        for line in f:
            pass
        last_workflow_id = line.rstrip()
    print(f"Last recorded wokflow ID was {last_workflow_id}")

    print(f"Creating new workflowIds.txt in {work_dir} directory using {cromshell_database}")
    workflow_database = pd.read_csv(cromshell_database, sep='\t')
    index_match = workflow_database[workflow_database["RUN_ID"] == last_workflow_id]["RUN_ID"].index[0]
    job_ids = workflow_database["RUN_ID"][ (1 + index_match): ].tolist()
    with open(workflow_id_tmp_file, "w") as output:
        output.writelines(id + '\n' for id in job_ids)

    print("Deleting previous workflowIds.txt file and replace with temp workflow_id")
    if os.path.exists(workflow_id_file):
        os.remove(workflow_id_file)
        copyfile(workflow_id_tmp_file, workflow_id_file)
        os.remove(workflow_id_tmp_file)
else:
    db = pd.read_csv(cromshell_database, sep='\t')
    job_ids = db["RUN_ID"].tolist()
    with open(workflow_id_file, "a") as output:
        output.writelines(id + '\n' for id in job_ids)

#########################################
# Run command to upload metadata to bq. #
#########################################

my_env = os.environ.copy()
my_env["CROMWELL_BASEURL"] = cromwell_base_url
my_env["GCP_PROJECT"] = gcp_project
my_env["DATASET_ID"] = "cromwell_monitoring"

metadata_upload_command = f"cromwell_metadata_bq < {workflow_id_file}"

print("\n## Upload the cromshell metadata ##")
print("Env. vars. set to:")
print(my_env)
print("Running the following command:")
print(metadata_upload_command)

subprocess.Popen(metadata_upload_command, cwd = work_dir, env = my_env)
