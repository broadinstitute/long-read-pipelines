We document here how runtime metrics and metadata are collected in a single BigQuery dataset (BQDS).

## Prerequisite and assumptions

We assume the the job was submitted using [cromshell](https://github.com/broadinstitute/cromshell) and controlled under a cromwell server.

--------------------------------------
## How does it work? Quick glance

### Collecting runtime metrics
The concurrent monitoring relies on a monitoring docker that can be specified via a workflow option when the job is submitted. The collected data (e.g. CPU/memory/disk percentage, etc) will be automatically sent to a BigQuery dataset. 

### Collecting job metadata
Metadata are relatively less volatile and are stored on the cromwell server that controlled the execution of the job. The metadata will be synced up to the same BQDS when the job is done.

### And ...
Much of these two functionalities have been implemented the repo [cromwell-task-monitor-bq](https://github.com/broadinstitute/cromwell-task-monitor-bq). We customize on top of that.

--------------------------------------
## Automated resources monitoring via WDL options

Users are expected to run an automated monitoring image alongside their job when it is launched using cromshell.

  * Specify `monitoring_image` in the workflow options json accompanying the submission (see [doc](https://cromwell.readthedocs.io/en/stable/wf_options/Google/)).

  * Launch the job via cromshell as you normally would
    
```bash
$ cromshell submit <WORKFLOW.wdl> <INPUT.json> <WORKFLOW_OPTIONS.json> <WORKFLOW.dependencies.zip>
```

When the monitoring image is used for the first time, it will automatically create a BQ dataset&mdash;named `cromwell_monitoring`&mdash;under the billing project the workflow is being run. The dataset has two date/timestamp partitioned tables&mdash;`metrics` and `runtime`.

  * `metrics` logs resources usage at 1 second interval (i.e. a time-series), whereas
  * `runtime` records the _static_ configurations of a particular VM instance for each task/shard (note that this is not necessarily the same as the `Runtime` specification that the WDL tasks requested, i.e. GCP does not necessarily provide exactly what the WDL asks for.)


--------------------------------------
## (Semi-)Automated Cromwell metadata upload
 
If the cromwell server that controls the job were configured appropriately, metadata should be sent to the BQDS `cromwell_monitoring` automatically once a job is finished.
However sometimes that server configuration is beyond our control, hence we need to upload the metadata from the cromwell server to the dataset, separately.

To setup, clone the [cromwell-task-monitor-bq](https://github.com/broadinstitute/cromwell-task-monitor-bq) repo. Then in its sub-directory "metadata/submit", compile a binary by `go build -o cromwell_metadata_bq`. Move the binary to a place in you `$PATH`.

A `metadata` table will appear under the BQDS `cromwell_monitoring` after either of the following commands ran for the first time. 


### manual upload
The metadata can be uploaded to the BQDS manually, given appropriate environment variables
 
```bash
$ export CROMWELL_BASEURL=<...>
$ export GCP_PROJECT=<...>
$ export DATASET_ID=cromwell_monitoring
```
Then

```bash
$ cromwell_metadata_bq < workflow_ids.txt
```

where "workflow_ids.txt" holds a list of workflow-ids (and nothing else) that have "succeeded" (according to `cromshell status`).

### automated upload using a python script
To automate the chore of checking execution status of jobs and pushing the metadata to the BQ dataset `cromwell_monitoring`, we offer a python3 script "cromwell.metadata.2BQ.py" (in the "utilities" sub-directory) that uploads the metadata of to a BQDS&mdash;preferably the same BQDS as above. We intend to have a cron job in the background that will periodically check the status of recently submitted jobs, and upload their metadata when possible.

First set the user-configurable variables in the script:

| variable    | value                                                                  |
|-------------|------------------------------------------------------------------------|
| gcp_project | the google cloud billing project                                       |
| work_dir    | working directory, assumed to hold the "cromwell-task-monitor-bq" repo |

Note that for the python script to work in cron, it needs to be run on the same machine where the job was launched with cromshell, because the script scans the jobs stored in 
"~/.cromshell/all.workflow.database.tsv" (which is automatically updated every time a job is submitted using cromshell) to search for finished jobs.
