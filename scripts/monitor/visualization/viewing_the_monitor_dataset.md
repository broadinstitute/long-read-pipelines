----------
## NOTE

When we initially experimented with visualizing BQ datasets by manipulating the BQ itself in the DB, and querying in one go, we ran into a probelm where the data is too big to be handled (for a workflow that had ~ 3K shards). We landed on a different way of querying the database&mdash;instead of creating a unified view as demonstrated below, query separately and batch query the huge `metrics` data.

Below is our experience in creating the unified view, which will be useful when the workflow isn't that humongous.


......

One needs to have an intuitive understandings of _dataset_, _tables_, _time-partitioned tables_, and _view_ of BigQuery before proceeding. It won't take much time.

......

Alright.

The repo [cromwell-task-monitor-bq](https://github.com/broadinstitute/cromwell-task-monitor-bq) forms the foundation of our current implementations, as we said in `collecting_metadata_and_runtime_metrics.md`. And we rely on the dataset `broad-dsde-methods.cromwell-monitoring`, composed of three tables `runtime`, `metrics` and `metadata`, generated in a way described in that document.


----------
## Convenience from an overall "view" of all three tables

To create a uniformed view of all three tables, we've generated a view named `monitor_everything_back_7days` using the following SQL query (note that the view is a logical view, not a materialized view, so every time a query is made to the view, the results are guaranteed to be the latest)

```sql
SELECT
  metadata.attempt AS meta_attempt,
  metadata.cpu_count AS meta_cpu,
  metadata.disk_mounts AS meta_disk_mounts,
  metadata.disk_total_gb AS meta_disk_total_gb,
  metadata.disk_types AS meta_disk_types,
  metadata.docker_image AS meta_docker_image,
  metadata.end_time AS meta_end_time,
  metadata.execution_status AS meta_execution_status,
  metadata.inputs AS meta_inputs,
  metadata.instance_name AS meta_instance_name,
  metadata.mem_total_gb AS meta_mem_total_gb,
  metadata.preemptible AS meta_preemptible,
  metadata.project_id AS meta_project_id,
  metadata.shard AS meta_shard,
  metadata.start_time AS meta_start_time,
  metadata.task_call_name AS meta_task_call_name,
  metadata.workflow_id AS meta_workflow_id,
  metadata.workflow_name AS meta_workflow_name,
  metadata.zone AS meta_zone,

  metrics.cpu_used_percent AS metrics_cpu_used_percent,
  metrics.disk_read_iops AS metrics_disk_read_iops,
  metrics.disk_used_gb AS metrics_disk_used_gb,
  metrics.disk_write_iops AS metrics_disk_write_iops,
  metrics.instance_id AS metrics_instance_id,
  metrics.mem_used_gb AS metrics_mem_used_gb,
  metrics.timestamp AS metrics_timestamp,

  runtime.attempt AS runtime_attempt,
  runtime.cpu_count AS runtime_cpu_count,
  runtime.cpu_platform AS runtime_cpu_platform,
  runtime.disk_mounts AS runtime_disk_mounts,
  runtime.disk_total_gb AS runtime_disk_total_gb,
  runtime.instance_id AS runtime_instance_id,
  runtime.instance_name AS runtime_instance_name,
  runtime.mem_total_gb AS runtime_mem_total_gb,
  runtime.preemptible AS runtime_preemptible,
  runtime.project_id AS runtime_project_id,
  runtime.shard AS runtime_shard,
  runtime.start_time AS runtime_start_time,
  runtime.task_call_name AS runtime_task_call_name,
  runtime.workflow_id AS runtime_workflow_id,
  runtime.zone AS runtime_zone

FROM
       `broad-dsde-methods.cromwell_monitoring.runtime`  runtime 
  JOIN `broad-dsde-methods.cromwell_monitoring.metrics`  metrics  ON runtime.instance_id = metrics.instance_id
  JOIN `broad-dsde-methods.cromwell_monitoring.metadata` metadata ON runtime.instance_name = metadata.instance_name

WHERE
      DATE(runtime.start_time) >= DATE_SUB(CURRENT_DATE(), INTERVAL 20 DAY)
  AND DATE(metrics.timestamp) >= DATE_SUB(CURRENT_DATE(), INTERVAL 20 DAY)
  AND DATE(metadata.start_time) >= DATE_SUB(CURRENT_DATE(), INTERVAL 20 DAY)

```

So yes, effectively we run a `SELECT *` into all three tables, but only going back 7 days (we think this is a good balance). Note that Google specifically said that using `LIMIT` does NOT save any money.

## Visualizing the data (3 sources)

We choose python to visualize the data from three sources, and compose a report.

### Retrieving the BQ table specific to a job

```python
import sys

import datetime
from datetime import timedelta

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import re

import numpy as np
import pandas as pd

from google.cloud import bigquery

sql_query = f"""
SELECT
  *
FROM
  broad-dsde-methods.cromwell_monitoring.monitor_everything_back_7days
WHERE
  runtime.workflow_id IN ({formated_workflow_ids})
"""

# define workflow and sub-workflow IDs
workflow_ids= ["\"6e0a7851-8618-4c40-a04d-a6c0ce125be6\"", 
               "\"2947c1f7-1a84-4724-abb8-c7c25c047dd5\"", 
               "\"93b8224b-e7f1-4226-adfa-01f3da15373d\""]
formated_workflow_ids = ','.join(workflow_ids)

# SQL into the view created above, restricting to the job we are interested in
client                   = bigquery.Client()
big_query_result         = client.query(sql_query)
runtime_metrics_metadata = big_query_result.to_dataframe()

```
