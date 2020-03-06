This python monitoring script will upload the cromwell metadata to BigQueary. The script uses the workflow IDs stored in the 
./cromshell tsv files to search the executed jobs in the cromwell server. It can be added as a cron job

In the python script edit the variables below:
**cromwellBaseUrl** = the Cromwell Url server being used
**gcpProject** = the google project being used when executing workflow
**cromwellTaskMonitorBqDirectory** = The location of the cromwell-task-monitor-bq directory that was cloned. specifically within cromwell-task-monitor-bq/metadata/submit 

## How does it work? 
A user will run a workflow using cromshell and the monitoring option will automatically run an image alongside each task and collect data. The collected data will be automatically sent to bigquery. On the same machine the user ran the cromshell, the user will have a cron job running that will periodically check the status of a submitted job, once that job stops then cron job will execute a command to obtain the runtime information for that specified workflow ID and push this to BQ.

## Workflow Monitoring
Users are expected to run an automated monitoring image along side their workflow using cromshell. 
1. Before running the workflow, in the options json add in the `monitoring_image` parameter for cromwell. For this image you will use the cromwell-bq docker image: `bshifaw/cromwell-task-monitor-bq:latest` (this should be moved to gcr for reliability soon)
2. Run cromshell as you normally would
`cromshell submit ValidateBam.wdl ValidateBam.input.json options/onlyMonitorImage.json ValidateBam.wdl.dependencies.zip`

 When initially running with this default image, it will automatically create tables in BigQueary within the project the workflow is being run. The table will be called `cromwell_monitoring`. Cromwell monitoring will have two partitioned tables `metrics` and `runtime`. The `metrics` table shows the resources usage at 1 sec interval time points, the `runtime` table records the VM instance resource . An additional table can be made by clicking on `Compose New Query` and entering in the following 
``` 
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
    `broad-dev-denis.cromwell_monitoring.metrics`
  WHERE
    DATE(timestamp) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
  GROUP BY
    instance_id
)SELECT
  r.project_id, r.zone, r.preemptible,
  r.workflow_id, r.task_call_name, r.shard, r.attempt,
  r.start_time runtime_start_time, runtime_duration_sec,
  cpu_platform, r.cpu_count, cpu_used_percent_avg,
  r.mem_total_gb, mem_used_gb_max,
  r.disk_mounts, r.disk_total_gb,
  (SELECT ARRAY_AGG(x IGNORE NULLS) FROM UNNEST(disk_used_gb_max) x) disk_used_gb_max,
  (SELECT ARRAY_AGG(x IGNORE NULLS) FROM UNNEST(disk_read_iops_avg) x) disk_read_iops_avg,
  (SELECT ARRAY_AGG(x IGNORE NULLS) FROM UNNEST(disk_write_iops_avg) x) disk_write_iops_avg
FROM
  `broad-dev-denis.cromwell_monitoring.runtime` r
JOIN
  metrics
USING (instance_id)
WHERE
  DATE(r.start_time) >= DATE_SUB(CURRENT_DATE(), INTERVAL 30 DAY)
ORDER BY
  r.start_time DESC
```
Then click on “Save View”, save it as you like i saved as runtime_metrics. This will create a view under the cromwell_monitoring table, the table combines the info from the metrics and runtime table create gather the useful info from both.
You can Query the view as such
```
SELECT *
FROM `broad-dsde-methods.cromwell_monitoring.runtime_metrics`
LIMIT 1000
```
## Cromwell MetaData using Cron Job

Initially you’ll need to clone the https://github.com/broadinstitute/cromwell-task-monitor-bq 
Change directories into to cromwell-task-monitor-bq/metadata/submit 
Run the following command 
`go build -o cromwell_metadata_bq`
Run this to upload the metadata to bigquery 
`CROMWELL_BASEURL=https://cromwell-v47.dsde-methods.broadinstitute.org GCP_PROJECT=<...> DATASET_ID=cromwell_monitoring ./cromwell_metadata_bq < workflow_ids.txt `
Where workflow_ids.txt is a list of workflow-ids in this form:
```
a5c1e1e7-76a5-4df0-9035-ef1996c17f0f
f7e339c9-1fa9-449d-8390-ee0cc4c39c5b
89fa4791-851a-4311-aa62-6864a2ef5fdf
```
After running for the first time the “metadata” dataset will appear under the cromwell_monitoring tab in BQ.

Next create a view using the following SQL command 
```
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

