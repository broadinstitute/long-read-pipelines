**NOTE!!!**

This instruction only works for **successfully** executed runs. The billing database is updated periodically (~daily). To guarantee accuracy, wait **24 hours** after the workflow
completes before running the query, as there is a risk you retrieve only partial results. Double check that costs for 
all expected tasks are present.

## How to get billed workflow costs:

1. Open BigQuery
2. Below the query editor box, click on "More" -> "Query setttings", then in the poped open side panel, choose "Legacy" under SQL dialect (because the SQL query shown here below is in legacy dialect), then "save"
3. Copy and paste this into the query box:

    ```sql
    SELECT
      GROUP_CONCAT(labels.key) WITHIN RECORD AS labels_key,
      GROUP_CONCAT(labels.value) WITHIN RECORD labels_value,
      sku.description,
      cost,
    FROM
      [broad-dsde-methods:Methods_billing_dump.gcp_billing_export_v1_009C7D_923007_219A6F]
    WHERE
      labels.key IN ("cromwell-workflow-id",
        "cromwell-workflow-name",
        "cromwell-sub-workflow-name",
        "wdl-task-name",
        "wdl-call-alias",
        "goog-dataproc-cluster-uuid",
        "goog-gke-node")
    HAVING
      labels_value LIKE "%f95259ff-ce19-4272-826f-cac18fa692b3%";
    ```

4. Replace the hash in the last line with your workflow ID
5. Hit RUN QUERY
6. Download as JSON
7. Source the script `func_big_query_cromwell_cost.R` and provide the appropriate arguments to the function to get summarized costs, broken down by task and SKU.