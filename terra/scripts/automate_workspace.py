#!/usr/bin/env python

import os
import re
import argparse

from dateutil.parser import parse

import pandas as pd
import firecloud.api as fapi
import json

from tinydb import TinyDB, Query


def main():
    parser = argparse.ArgumentParser(
        description='Automatically run workflows when new data is added to a workspace',
        prog='automate_workspace'
    )
    parser.add_argument('-n', '--namespace', required=True, type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', required=True, type=str, help="Terra workspace")
    parser.add_argument('-m', '--min-date', default='2021-05-01T00:00:00', type=str, help="Lower date cutoff")
    args = parser.parse_args()

    allowed_states = ['Failed', 'Aborted']
    min_date = parse(args.min_date)

    ent = fapi.get_entities(args.namespace, args.workspace, 'sample').json()
    if len(ent) > 0:
        tbl = pd.DataFrame(list(map(lambda e: e['attributes'], ent)))
        tbl["entity:sample_id"] = list(map(lambda f: f['name'], ent))

        # current_config = fapi.get_workspace_config(args.namespace, args.workspace, args.namespace, 'PBFlowcell')
        # current_config_json = current_config.json()
        # print(json.dumps(current_config_json, indent=4, sort_keys=True))

        s = {}
        submissions = fapi.list_submissions(args.namespace, args.workspace).json()
        for submission in submissions:
            workflowStatus = list(submission['workflowStatuses'].keys())[0]
            sample_id = re.sub("_.+", "", submission['submissionEntity']['entityName'])

            if workflowStatus not in allowed_states:
                s[sample_id] = s.get(sample_id, 0) + 1

        for i in range(len(tbl)):
            flowcell_date = parse(re.sub("[TZ\.].*", "", tbl['created_at'][i]))
            sample_id = tbl["entity:sample_id"][i]

            if flowcell_date > min_date and sample_id not in s:
                response = fapi.create_submission(
                    args.namespace,
                    args.workspace,
                    args.namespace,
                    'PBFlowcell',
                    tbl["entity:sample_id"][i],
                    'sample',
                    use_callcache=True
                )
                submission_result = response.json()
                submission_id = submission_result['submissionId']

                print(f'sample_id={tbl["entity:sample_id"][i]} bio_sample={tbl["bio_sample"][i]} well_sample={tbl["well_sample"][i]} experiment_type={tbl["experiment_type"][i]} submission_id={submission_id}')


if __name__ == "__main__":
    main()
