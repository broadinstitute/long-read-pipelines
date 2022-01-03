#!/usr/bin/env python

import os
import re
import argparse

from dateutil.parser import parse
from numpy import exp

import pandas as pd
import firecloud.api as fapi
import json


def main():
    parser = argparse.ArgumentParser(
        description='Automatically run workflows when new data is added to a workspace',
        prog='automate_workspace'
    )
    parser.add_argument('-n', '--namespace', required=True, type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', required=True, type=str, help="Terra workspace")
    parser.add_argument('-m', '--min-date', default='2021-12-12T00:00:00Z', type=str, help="Lower date cutoff")
    parser.add_argument('-b', '--branch', default=['main'], action='extend', nargs="+", help="The branch to require for automated job submissions")
    parser.add_argument('-r', '--run', action='store_true', help="Turn off the default dry-run mode")
    args = parser.parse_args()

    # Get a table of previous or current job submissions
    job_submissions = get_workflow_submission_statuses(args.namespace, args.workspace)

    # Get list of flowcells (we'll need this to determine the data type at the sample_set level)
    s_tbl = load_table(args, 'sample')
    s_tbl['created_at']= pd.to_datetime(s_tbl['created_at'])

    workflows = {
        'PBFlowcell': set(['CCS', 'CLR', 'ISOSEQ']),
        'PBCCSWholeGenome': set(['CCS']),
        'PBCLRWholeGenome': set(['CLR']),
        'PBCCSIsoSeq': set(['ISOSEQ']),
    }

    # Iterate through workflows looking for jobs to run
    for wf in workflows:
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
                    tbl_w_jobs_subset = tbl_w_jobs[tbl_w_jobs[f'entity:{root_entity_type}_id'] == id].sort_values('submission_date', ascending=False)

                    samples = None
                    created_at = None
                    experiment_type = None
                    prereq_jobs = []
                    if len(tbl_w_jobs_subset) > 0:
                        if 'samples' in tbl_w_jobs_subset and len(tbl_w_jobs_subset['samples'][0]['items']) > 0:
                            s_tbl_subset = s_tbl[s_tbl['entity:sample_id'] == tbl_w_jobs_subset['samples'][0]['items'][0]['entityName']].sort_values('created_at', ascending=False)
                            samples = list(s_tbl_subset['entity:sample_id'])
                            created_at = list(s_tbl_subset['created_at'])[0]
                            experiment_type = list(s_tbl_subset['experiment_type'])[0]

                            prereq_jobs = job_submissions[(job_submissions['sample_id'].isin(samples)) &
                                                          (job_submissions['method_configuration_name'] == 'PBFlowcell') &
                                                          (job_submissions['workflow_status'] == 'Succeeded')].sort_values('submission_date', ascending=False)
                        elif 'created_at' in tbl_w_jobs_subset and 'experiment_type' in tbl_w_jobs_subset:
                            samples = list(tbl_w_jobs_subset['entity:sample_id'])
                            created_at = list(tbl_w_jobs_subset['created_at'])[0]
                            experiment_type = list(tbl_w_jobs_subset['experiment_type'])[0]
                            prereq_jobs = [None]

                    process = False
                    if experiment_type in workflows[wf]:
                        if len(tbl_w_jobs_subset) == 0:
                            process = True
                        else:
                            has_succeeded = tbl_w_jobs_subset['workflow_status'][0] == 'Succeeded'
                            is_running = 'Running' in list(tbl_w_jobs_subset['workflow_status'])
                            has_exceeded_failure_limit = len(list(filter(lambda x: x != "Succeeded" and x != "Running", list(tbl_w_jobs_subset['workflow_status'])))) >= 3
                            created_before_min_date = pd.to_datetime(created_at) < pd.to_datetime(args.min_date)
                            prereqs_satisfied = len(prereq_jobs) > 0

                            process = not has_succeeded and not is_running and not has_exceeded_failure_limit and not created_before_min_date and prereqs_satisfied

                    if process:
                        print(f'{wf} {id} {process} {created_at}')

                        submission_id = "dry-run"
                        submission_result = ""
                        entity_id = list(set(tbl_w_jobs_subset[f'entity:{root_entity_type}_id']))[0]
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



    # current_config = fapi.get_workspace_config(args.namespace, args.workspace, args.namespace, 'PBFlowcell')
    # current_config_json = current_config.json()
    # print(json.dumps(current_config_json, indent=4, sort_keys=True))

    # samples_to_process = get_samples_to_process(args)
    # min_date = parse(args.min_date)

    # allowed_states = ['Failed', 'Aborted']

    # # Load all available samples from workspace
    # ent = fapi.get_entities(args.namespace, args.workspace, 'sample').json()
    # if len(ent) > 0:
    #     tbl = pd.DataFrame(list(map(lambda e: e['attributes'], ent)))
    #     tbl["entity:sample_id"] = list(map(lambda f: f['name'], ent))

    #     # Get a list of previous or current job submissions
    #     processed_samples = {}
    #     submissions = fapi.list_submissions(args.namespace, args.workspace).json()
    #     for submission in submissions:
    #         workflowStatus = list(submission['workflowStatuses'].keys())[0]
    #         sample_id = re.sub("_.+", "", submission['submissionEntity']['entityName'])

    #         if workflowStatus not in allowed_states:
    #             processed_samples[sample_id] = processed_samples.get(sample_id, 0) + 1

    #     for i in range(len(tbl)):
    #         flowcell_date = parse(re.sub("[TZ\.].*", "", tbl['created_at'][i]))
    #         sample_id = tbl["entity:sample_id"][i]

    #         if flowcell_date > min_date or sample_id in samples_to_process:
    #             if sample_id in processed_samples and sample_id not in samples_to_process:
    #                 if args.verbose:
    #                     print(f'sample_id={tbl["entity:sample_id"][i]} bio_sample={tbl["bio_sample"][i]} well_sample={tbl["well_sample"][i]} experiment_type={tbl["experiment_type"][i]} submission_id=done')
    #             else:
    #                 submission_id = "dry-run"
    #                 if args.run:
    #                     response = fapi.create_submission(
    #                         args.namespace,
    #                         args.workspace,
    #                         args.namespace,
    #                         'PBFlowcell',
    #                         tbl["entity:sample_id"][i],
    #                         'sample',
    #                         use_callcache=True
    #                     )
    #                     submission_result = response.json()
    #                     submission_id = submission_result['submissionId']

    #                 print(f'sample_id={tbl["entity:sample_id"][i]} bio_sample={tbl["bio_sample"][i]} well_sample={tbl["well_sample"][i]} experiment_type={tbl["experiment_type"][i]} submission_id={submission_id}')


def get_samples_to_process(args):
    samples_to_process = set()
    if args.sample_ids is not None:
        for s in args.sample_ids:
            if os.path.isfile(s):
                with open(s) as fp:
                    lines = fp.readlines()
                    for line in lines:
                        samples_to_process.add(line.strip())
            else:
                samples_to_process.add(s)
    return samples_to_process


def get_workflow_submission_statuses(namespace, workspace):
    tbl_header = ["method_configuration_name", "sample_id", "workflow_status", "submission_date"]
    tbl_rows = []
    submissions = fapi.list_submissions(namespace, workspace).json()
    for submission in submissions:
        method_configuration_name = re.sub("_(?:.(?!_))+$", "", submission['methodConfigurationName'])
        sample_id = re.sub("_(?:.(?!_))+$", "", submission['submissionEntity']['entityName'])
        submission_date = submission['submissionDate']

        for workflow_status in filter(lambda x: x != "Aborted", submission['workflowStatuses']):
            for _ in range(0, int(submission['workflowStatuses'][workflow_status])):
                tbl_rows.append([method_configuration_name, sample_id, workflow_status, submission_date])
    
    tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)
    tbl_new['submission_date']= pd.to_datetime(tbl_new['submission_date'])

    return tbl_new


if __name__ == "__main__":
    main()
