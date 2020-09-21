We document here how runtime metrics and metadata are collected in a single BigQuery dataset (BQDS).

## Prerequisite and assumptions

We assume the the job was submitted using [cromshell](https://github.com/broadinstitute/cromshell) and controlled under a cromwell server.

--------------------------------------
## How does it work? Quick glance

### Collecting runtime metrics
The concurrent monitoring relies on a monitoring docker that can be specified via a workflow option when the job is submitted. The collected data (e.g. CPU/memory/disk percentage, etc) will be automatically sent to a BigQuery dataset. 

### Collecting job metadata
Metadata are relatively less volatile and are stored in the cromwell server that controlled the execution of the job. The metadata will be synced up to the same BQDS when the job is done.


Much of these two functionalities have been implemented in [this repo](https://github.com/broadinstitute/cromwell-task-monitor-bq). We customize on top of that.

--------------------------------------
## Workflow resources monitoring

Users are expected to run an automated monitoring image along side their job when it is launched using cromshell.

  1. Before running the workflow, specify `monitoring_image` in the workflow options json accompanying the submission (see [doc](https://cromwell.readthedocs.io/en/stable/wf_options/Google/)). An image known to be working is `us.gcr.io/broad-dsde-methods/cromwell-task-monitor-bq:latest` (this should be moved to gcr for reliability soon). 

  2. Launch the job via cromshell as you normally would
    
```bash
$ cromshell submit <WORKFLOW.wdl> <INPUT.json> <WORKFLOW_OPTIONS.json> <WORKFLOW.dependencies.zip>
```

  3. When the monitoring image is used for the first time, it will automatically create a BQ dataset&mdash;named `cromwell_monitoring`&mdash;under the billing project the workflow is being run. The dataset has two partitioned tables `metrics` and `runtime`. `metrics` logs resources usage at 1 second interval (time-series), whereas `runtime` records the _static_ specifications of a particular VM instance for each task/shard. 

  4. To view the two tables together, an additional table can be made by clicking on "Compose New Query" and entering the following 
 
```sql
WITH
  metrics AS (
  SELECT
    instance_id,
    TIMESTAMP_DIFF(MAX(timestamp), MIN(timestamp), SECOND) runtime_duration_sec,
    AVG((
      SELECT
        AVG(p)
      FROM
        UNNEST(cpu_used_percent) p)) cpu_used_percent_avg,
    MAX(mem_used_gb) mem_used_gb_max,
    [MAX(disk_used_gb[
    OFFSET
      (0)]),
    MAX(disk_used_gb[SAFE_OFFSET(1)])] disk_used_gb_max,
    [AVG(disk_read_iops[
    OFFSET
      (0)]),
    AVG(disk_read_iops[SAFE_OFFSET(1)])] disk_read_iops_avg,
    [AVG(disk_write_iops[
    OFFSET
      (0)]),
    AVG(disk_write_iops[SAFE_OFFSET(1)])] disk_write_iops_avg
  FROM
    `broad-dev-denis.cromwell_monitoring.metrics`
  WHERE
    DATE(timestamp) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
  GROUP BY
    instance_id )
SELECT
  r.project_id,
  r.zone,
  r.preemptible,
  r.workflow_id,
  r.task_call_name,
  r.shard,
  r.attempt,
  r.start_time runtime_start_time,
  runtime_duration_sec,
  cpu_platform,
  r.cpu_count,
  cpu_used_percent_avg,
  r.mem_total_gb,
  mem_used_gb_max,
  r.disk_mounts,
  r.disk_total_gb,
  (
  SELECT
    ARRAY_AGG(x IGNORE NULLS)
  FROM
    UNNEST(disk_used_gb_max) x) disk_used_gb_max,
  (
  SELECT
    ARRAY_AGG(x IGNORE NULLS)
  FROM
    UNNEST(disk_read_iops_avg) x) disk_read_iops_avg,
  (
  SELECT
    ARRAY_AGG(x IGNORE NULLS)
  FROM
    UNNEST(disk_write_iops_avg) x) disk_write_iops_avg
FROM
  `broad-dev-denis.cromwell_monitoring.runtime` r
JOIN
  metrics
USING
  (instance_id)
WHERE
  DATE(r.start_time) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
ORDER BY
  r.start_time DESC
```
   Click on “Save View”, save it as you like (e.g. "runtime_metrics'). This will create a view under the `cromwell_monitoring` dataset&mdash;the table combines the info from the `metrics` and `runtime` table and allows you gather the useful info from both. You can then query the view as such

```sql
    SELECT *
    FROM `broad-dsde-methods.cromwell_monitoring.runtime_metrics`
    LIMIT 1000
```

--------------------------------------
## Cromwell metadata upload
 
If the cromwell server that controls the job were configured appropriately, metadata should be sent to the BQDS `cromwell_monitoring` automatically once a job is finished.
However sometimes that server configuration is beyond our control, hence we need to upload the metadata from the cromwell server to the dataset, separately.

To setup, clone the [cromwell-task-monitor-bq](https://github.com/broadinstitute/cromwell-task-monitor-bq) repo. Then in its sub-directory "metadata/submit", compile a binary by `go build -o cromwell_metadata_bq`. Move the binary to a place in you `$PATH`.

### manual upload
The metadata can be uploaded to the BQDS manually via providing the appropriate environment variables
 
```bash
$ CROMWELL_BASEURL=<...>
$ GCP_PROJECT=<...>
$ DATASET_ID=cromwell_monitoring
$ cromwell_metadata_bq < workflow_ids.txt
```
where "workflow_ids.txt" holds a list of workflow-ids (and nothing else) that have "succeeded" (according to `cromshell status`).

### automated upload using a python script
To automate the chore of checking execution status of jobs and pushing the metadata to the BQ dataset `cromwell_monitoring`, we offer a python3 script "cromwell.metadata.2BQ.py" that uploads the metadata of to a BQDS&mdash;preferably the same BQDS as above. We intend to have a cron job in the background that will periodically check the status of recently submitted jobs, and upload their metadata when possible.

First set the user-configurable variables in the script:

| variable    | value                                                                  |
|-------------|------------------------------------------------------------------------|
| gcp_project | the google cloud billing project                                       |
| work_dir    | working directory, assumed to hold the "cromwell-task-monitor-bq" repo |

Note that for the python script to work in cron, it needs to be run on the same machine where the job was launched with cromshell, because the script scans the jobs stored in 
"~/.cromshell/all.workflow.database.tsv" (which is automatically updated every time a job is submitted using cromshell) to search for finished jobs.

### querying the "metadata" table

A “metadata” table will appear under the BQDS `cromwell_monitoring` after the above command ran for the first time. 

Next create a view using the following SQL command 

```sql
WITH metrics AS (
  SELECT
	instance_id,
	TIMESTAMP_DIFF(MAX(timestamp), MIN(timestamp), SECOND) runtime_duration_sec,
	AVG((SELECT AVG(p) FROM UNNEST(cpu_used_percent) p)) cpu_used_percent_avg,
	MAX(mem_used_gb) mem_used_gb_max,
	[MAX(disk_used_gb[OFFSET(0)]), MAX(disk_used_gb[SAFE_OFFSET(1)])] disk_used_gb_max,
	[AVG(disk_read_iops[OFFSET(0)]), AVG(disk_read_iops[SAFE_OFFSET(1)])] disk_read_iops_avg,
	[AVG(disk_write_iops[OFFSET(0)]), AVG(disk_write_iops[SAFE_OFFSET(1)])] disk_write_iops_avg
  FROM
	`broad-epi-dev.cromwell_monitoring.metrics`
  WHERE
	DATE(timestamp) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
  GROUP BY
	instance_id
)SELECT
  r.project_id, r.zone, r.preemptible,
  r.workflow_id, workflow_name, r.task_call_name, r.shard, r.attempt, execution_status,
  m.start_time metadata_start_time, TIMESTAMP_DIFF(m.end_time, m.start_time, SECOND) metadata_duration_sec, runtime_duration_sec,
  cpu_platform, r.cpu_count, cpu_used_percent_avg,
  r.mem_total_gb, mem_used_gb_max,
  r.disk_mounts, disk_types, r.disk_total_gb,
  (SELECT ARRAY_AGG(x IGNORE NULLS) FROM UNNEST(disk_used_gb_max) x) disk_used_gb_max,
  (SELECT ARRAY_AGG(x IGNORE NULLS) FROM UNNEST(disk_read_iops_avg) x) disk_read_iops_avg,
  (SELECT ARRAY_AGG(x IGNORE NULLS) FROM UNNEST(disk_write_iops_avg) x) disk_write_iops_avg,
  docker_image, inputs
FROM
  `broad-epi-dev.cromwell_monitoring.runtime` r
JOIN
  metrics
USING (instance_id)
JOIN
  `broad-epi-dev.cromwell_monitoring.metadata` m
USING (instance_name)
WHERE
  DATE(r.start_time) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
  AND
  DATE(m.start_time) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
ORDER BY
  r.start_time DESC
```

