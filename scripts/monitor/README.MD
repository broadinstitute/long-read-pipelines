This python script will upload the cromwell metadata to BigQueary. It can be added as a cron job

In the python script edit the variables below

#cromwell-task-monitor-bq Variables
cromwellBaseUrl = the Cromwell Url server being used
gcpProject = the google project being used when executing workflow
cromwellTaskMonitorBqDirectory = The location of the cromwell-task-monitor-bq directory that was cloned. specifically within cromwell-task-monitor-bq/metadata/submit 

