#!/usr/bin/env python

import os
import re
import argparse

import pandas as pd
import firecloud.api as fapi
import json

from tinydb import TinyDB, Query


def main():
    parser = argparse.ArgumentParser(
        description='Automatically run workflows when new data is added to a workspace',
        prog='automate_workspace'
    )
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    args = parser.parse_args()

    db = TinyDB(f'{args.namespace}.{args.workspace}.json')

    current_config = fapi.get_workspace_config(args.namespace, args.workspace, args.namespace, 'PBFlowcell')
    current_config_json = current_config.json()

    print(json.dumps(current_config_json, indent=4, sort_keys=True))

    ent = fapi.get_entities(args.namespace, args.workspace, 'sample').json()
    if len(ent) > 0:
        tbl = pd.DataFrame(list(map(lambda e: e['attributes'], ent)))
        tbl["entity:sample_id"] = list(map(lambda f: f['name'], ent))

        print(tbl)

        # for i in range(len(tbl)):
        #     print(tbl.iloc(i))

    # print(json.dumps(fapi.get_submission_queue().json(), indent=4, sort_keys=True))
    # print(json.dumps(fapi.list_submissions(args.namespace, args.workspace).json(), indent=4, sort_keys=True))

    submissions = fapi.list_submissions(args.namespace, args.workspace).json()
    for s in submissions:
        print(f"{s['submissionEntity']['entityName']} {s['workflowStatuses']}")


if __name__ == "__main__":
    main()


# {
#     "deleteIntermediateOutputFiles": true,
#     "methodConfigurationDeleted": true,
#     "methodConfigurationName": "PBFlowcell_0Vr3Q53xqoA",
#     "methodConfigurationNamespace": "production-long-reads",
#     "status": "Submitted",
#     "submissionDate": "2021-05-08T07:13:10.290Z",
#     "submissionEntity": {
#         "entityName": "01300ed1-77fe-41bb-92c2-7907dfb910c5",
#         "entityType": "sample"
#     },
#     "submissionId": "42a145b8-a335-4fa8-8d00-7e0dc1e3206c",
#     "submitter": "kiran@broadinstitute.org",
#     "useCallCache": true,
#     "workflowStatuses": {
#         "Running": 1
#     }
# }