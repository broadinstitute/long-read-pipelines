We aim at monitoring cromwell-controlled jobs using three different sources of data:

* the metadata that contain the VM specification and inputs of the job
* the runtime (i.e. time series) resource usage of each sub-workflow and tasks of the job
* the cost incurred

We collect these data from:

* the cromwell server that controls the execution of the job
* an output of a monitoring docker that was specified as an option when the job was launched
* a BigQuery database where the cost information is updated periodically

Currently the first two tasks are documented in "collecting\_metadata\_and\_runtime\_metrics.md", while cost data retrieval is documented in "collecting\_cost\_from\_big\_query.md".

(Note: We document viewing&mdash;in the SQL and visualization terms&mdash;the non-cost part of data in "viewing\_the\_monitor_dataset.md" when the workflow isn't too large to be handled by the API call to BQ.) 
