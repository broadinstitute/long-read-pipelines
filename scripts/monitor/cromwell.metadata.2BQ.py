#!/usr/bin/env python3
import os
from os import path
from shutil import copyfile
import csv
import time
import sys
import subprocess
import pandas as pd

#cromwell-task-monitor-bq Variables
cromwellBaseUrl = "https://cromwell-v47.dsde-methods.broadinstitute.org"
gcpProject = "broad-dsde-methods"
cromwellTaskMonitorBqDirectory = "/Users/bshifaw/work/Broadinstitute/cromwell-task-monitor-bq/metadata/submit/"




#obtain workflow id
print("/n" + "## Obtaining Wokflow IDs ##")
HOME = os.getenv('HOME')
cromshellDatabase = HOME + '/.cromshell/all.workflow.database.tsv'
pythonScriptPath = sys.path[0]
workflowID = pythonScriptPath + '/workflowIds.txt'
workflowIDTemp = pythonScriptPath + '/workflowIds_temp.txt'


######### read timefile
# what was the last time that was created 
# grep that time in the databasetsve and copy workflowids after that match. 
# if nothing was updated then exit with message.


if path.exists(workflowID) and os.stat(workflowID).st_size > 0:
	print("Grep last wokflow ID") 
	with open(workflowID, "r") as f:
	    for line in f:
	        pass
	    last_workflowID = line.rstrip()
	f.close
	print("Last recorded wokflow ID was " + last_workflowID)

	print("Creating new workflowIds.txt in " + pythonScriptPath + " directory using " + cromshellDatabase)
	
	workflow_database = pd.read_csv(cromshellDatabase, sep='\t')
    
	print(workflow_database.iloc[:,2].str.match(last_workflowID))

	index_match = workflow_database.index[workflow_database.iloc[:,2].str.match(last_workflowID)].tolist()[0]
	
	with open(workflowIDTemp, "w") as the_file:
		for id in workflow_database.iloc[(1+index_match):,2]:
			the_file.write(id+'\n')

	print("Deleting previous workflowIds.txt file and replace with temp workflowID")
	if path.exists(workflowID):
		os.remove(workflowID)
		copyfile(workflowIDTemp, workflowID)
		os.remove(workflowIDTemp)
else: 
	f = open(workflowID, "a")
	with open(cromshellDatabase) as tsvfile:
		tsvreader = csv.reader(tsvfile, delimiter="\t")
	
		for line in tsvreader:
			if line[2] != "RUN_ID":
				f.write(line[2] + "\n")
f.close

############ overwrite timefile to say this was last time i checked

#Run command to upload metadata to bq.
metadataUploadCommand = 'CROMWELL_BASEURL=' + cromwellBaseUrl + ' GCP_PROJECT=' + gcpProject + ' DATASET_ID=cromwell_monitoring ./cromwell_metadata_bq < ' + workflowID 
print("/n" + "## Upload the cromshell metadata ##")
print("Running the following command:")
print(metadataUploadCommand)
#shell=True is said to be a security risk, may need to change https://stackoverflow.com/questions/18962785/oserror-errno-2-no-such-file-or-directory-while-using-python-subprocess-in-dj
#subprocess.Popen(metadataUploadCommand, cwd=cromwellTaskMonitorBqDirectory, shell=True)




# useful code for timeing line[1:]-time.strftime('%Y%m%d_%H%M%S')
