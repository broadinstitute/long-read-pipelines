To view the two tables `runtime` and `metrics` under the BQ dataset `cromwell_monitoring` together, an additional table can be made by clicking on "Compose New Query" and entering the following 
 
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


To merge three tables together

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