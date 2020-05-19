**NOTE!!!**

This instruction only works for **successfully** executed runs. The billing database is updated periodically (~daily). To guarantee accuracy, wait **24 hours** after the workflow
completes before running the query, as there is a risk you retrieve only partial results. Double check that costs for 
all expected tasks are present.

## How to get billed workflow costs:

1. Open BigQuery
2. Copy and paste this into the query box:

    ```sql
     SELECT
      (SELECT value FROM UNNEST(labels) WHERE key = 'cromwell-workflow-id') AS workflow_id,
      (SELECT value FROM UNNEST(labels) WHERE key = 'cromwell-workflow-name') AS workflow_name,
      (SELECT value FROM UNNEST(labels) WHERE key = 'cromwell-sub-workflow-name') AS sub_workflow_id,
      (SELECT value FROM UNNEST(labels) WHERE key = 'wdl-task-name') AS task_name,
      (SELECT value FROM UNNEST(labels) WHERE key = 'wdl-call-alias') AS task_alias,
      (SELECT value FROM UNNEST(labels) WHERE key = 'goog-gke-node') AS node,
     
      sku.description AS sku_description,
      cost,
      cost_type,
      
      usage_start_time,
      usage_end_time,
      (SELECT value FROM UNNEST(system_labels) WHERE key = 'compute.googleapis.com/cores') AS cores,
      (SELECT value FROM UNNEST(system_labels) WHERE key = 'compute.googleapis.com/memory') AS memory,
      usage
    
    
    FROM
      `broad-dsde-methods.Methods_billing_dump.gcp_billing_export_v1_009C7D_923007_219A6F`
    
    LEFT JOIN
      UNNEST(labels) AS label
    
    WHERE
      cost > 0.0
      AND usage_start_time>='2020-05-08 00:00:00'
      AND usage_start_time<='2020-05-09 00:00:00'
      AND label.key IN ("cromwell-workflow-id",
                        "cromwell-workflow-name",
                        "cromwell-sub-workflow-name",
                        "wdl-task-name",
                        "wdl-call-alias",
                        "goog-gke-node")
      AND label.value LIKE "%93b8224b-e7f1-4226-adfa-01f3da15373d%"
    ```

3. Replace the start/end time with appropriate time stamp (leave some slack)
4. Replace the workflow ID in the last line with your workflow ID (DO NOT MOVE THIS LINE)
5. Hit RUN QUERY
6. Download as JSON (if your workflow is unrealistically big, the JSON file may be limited in row count)
7. Source the script `func_big_query_cromwell_cost.R` and provide the appropriate arguments to the function to get summarized costs, broken down by task and SKU.