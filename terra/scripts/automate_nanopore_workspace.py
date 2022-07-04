#!/usr/bin/env python

import os
import re
import argparse

from dateutil.parser import parse
from numpy import exp

import pandas as pd
import firecloud.api as fapi
import json


pd.set_option('max_columns', 2000)
pd.set_option('max_rows', 200)
pd.set_option("max_colwidth", None)


def main():
    parser = argparse.ArgumentParser(
        description='Automatically run workflows when new data is added to a Nanopore workspace',
        prog='automate_nanopore_workspace'
    )
    parser.add_argument('-n', '--namespace', required=True, type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', required=True, type=str, help="Terra workspace")
    parser.add_argument('-b', '--branch', default=['main'], action='extend', nargs="+", help="The branch to require for automated job submissions")
    parser.add_argument('-r', '--run', action='store_true', help="Turn off the default dry-run mode")
    args = parser.parse_args()

    # Get a table of previous or current job submissions
    job_submissions = get_workflow_submission_statuses(args.namespace, args.workspace)
    #print(job_submissions)

    # Get list of flowcells (we'll need this to determine the data type at the sample_set level)
    s_tbl = load_table(args, 'sample')
    try:
        s_tbl['created_at']= pd.to_datetime(s_tbl['created_at'])
    except:
        s_tbl['created_at']= pd.to_datetime('2019-01-01T00:00:00.000Z')

    workflow_dependencies = {
        #'ONTFlowcell': None,
        'ONTWholeGenome': ['ONTFlowcell'],
    }

    # Iterate through workflows looking for jobs to run
    for wf in workflow_dependencies:
        current_config = fapi.get_workspace_config(args.namespace, args.workspace, args.namespace, wf)
        if current_config.status_code != 404:
            current_config_json = current_config.json()

            method_version = current_config_json['methodRepoMethod']['methodVersion']
            if method_version in args.branch:
                root_entity_type = current_config_json['rootEntityType']

                tbl = load_table(args, root_entity_type)
                tbl_w_jobs = pd.merge(tbl.set_index(f'entity:{root_entity_type}_id', drop=False),
                                      job_submissions.set_index(f'sample_id', drop=True),
                                      left_index=True, right_index=True, how="outer")

                for id in tbl_w_jobs[f'entity:{root_entity_type}_id'].unique():
                    tbl_w_jobs_subset = tbl_w_jobs[(tbl_w_jobs[f'entity:{root_entity_type}_id'] == id) &
                                                   (tbl_w_jobs[f'method_configuration_name'] == wf)].sort_values('submission_date', ascending=False)

                    dependencies = workflow_dependencies[wf]

                    dependencies_satisfied = False
                    if dependencies is None:
                        dependencies_satisfied = True
                    else:
                        print(id)

                        for dependency in workflow_dependencies[wf]:
                            print(dependency)
                            print(job_submissions[(job_submissions['sample_id'] == id) & (job_submissions['method_configuration_name'] == dependency)])

                            prereq_jobs = job_submissions[(job_submissions['sample_id'] == id) &
                                                          (job_submissions['method_configuration_name'] == dependency) &
                                                          (job_submissions['workflow_status'] == 'Succeeded')].sort_values('submission_date', ascending=False)
                            dependencies_satisfied |= len(prereq_jobs) > 0

                    has_succeeded = False if len(tbl_w_jobs_subset) == 0 else tbl_w_jobs_subset['workflow_status'][0] == 'Succeeded'
                    is_running = 'Running' in list(tbl_w_jobs_subset['workflow_status'])
                    has_exceeded_failure_limit = len(list(filter(lambda x: x != "Succeeded" and x != "Running", list(tbl_w_jobs_subset['workflow_status'])))) >= 3

                    process = not has_succeeded and not is_running and not has_exceeded_failure_limit and dependencies_satisfied 

                    if process:
                        submission_id = "dry-run"
                        submission_result = ""
                        entity_id = id
                        if args.run:
                            response = fapi.create_submission(
                                args.namespace,
                                args.workspace,
                                args.namespace,
                                wf,
                                entity_id,
                                root_entity_type,
                                use_callcache=True
                            )
                            submission_result = response.json()
                            submission_id = submission_result['submissionId']
                            
                        print(f'namespace={args.namespace}, workspace={args.workspace}, wf={wf}, entity={entity_id}, submission_id={submission_id}')
            else:
                print(f'Skipping {wf} because current configuration is set to "{method_version}" rather than "{args.branch}".')


def load_table(args, root_entity_type):
    ent = fapi.get_entities(args.namespace, args.workspace, root_entity_type).json()
    tbl = pd.DataFrame(list(map(lambda e: e['attributes'], ent)))
    tbl[f"entity:{root_entity_type}_id"] = list(map(lambda f: f['name'], ent))

    return tbl


def get_workflow_submission_statuses(namespace, workspace):
    tbl_header = ["method_configuration_name", "sample_id", "workflow_status", "submission_date"]
    tbl_rows = []
    submissions = fapi.list_submissions(namespace, workspace).json()
    for submission in submissions:
        method_configuration_name = re.sub("_(?:.(?!_))+$", "", submission['methodConfigurationName'])
        # sample_id = re.sub("_(?:.(?!_))+$", "", submission['submissionEntity']['entityName'])
        sample_id = submission['submissionEntity']['entityName']
        submission_date = submission['submissionDate']

        for workflow_status in filter(lambda x: x != "Aborted", submission['workflowStatuses']):
            for _ in range(0, int(submission['workflowStatuses'][workflow_status])):
                tbl_rows.append([method_configuration_name, sample_id, workflow_status, submission_date])
    
    tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)
    tbl_new['submission_date']= pd.to_datetime(tbl_new['submission_date'])

    return tbl_new


if __name__ == "__main__":
    main()
