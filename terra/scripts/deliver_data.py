#!/usr/bin/env python

import os
import re
import math
import hashlib
import argparse
import subprocess

import numpy as np
import pandas as pd
import firecloud.api as fapi

from google.cloud import bigquery
from google.cloud import storage
from google.api_core.exceptions import NotFound

from collections import OrderedDict

import xmltodict
import pprint


pd.set_option('max_columns', 200)
pd.set_option('max_rows', 200)
pd.set_option("max_colwidth", None)


def copy_file(src, dst):
    # (bucket_name, blob_name) = re.split("/", re.sub("gs://", "", src), 1)
    # (destination_bucket_name, destination_blob_name) = re.split("/", re.sub("gs://", "", dst), 1)
    #
    # storage_client = storage.Client()
    #
    # source_bucket = storage_client.bucket(bucket_name)
    # source_blob = source_bucket.blob(blob_name)
    # destination_bucket = storage_client.bucket(destination_bucket_name)
    #
    # blob_copy = source_bucket.copy_blob(
    #     source_blob, destination_bucket, destination_blob_name
    # )
    # dst_blob.rewrite(src_blob)
    #
    # print(
    #     "Blob {} in bucket {} copied to blob {} in bucket {}.".format(
    #         source_blob.name,
    #         source_bucket.name,
    #         blob_copy.name,
    #         destination_bucket.name,
    #     )
    # )

    subprocess.run(["gsutil", "cp", "-n", src, dst])


def load_table(namespace, workspace, table_name, store_membership=False):
    ent_old = fapi.get_entities(namespace, workspace, table_name).json()
    tbl_old = None

    membership = None
    if len(ent_old) > 0:
        tbl_old = pd.DataFrame(list(map(lambda e: e['attributes'], ent_old)))
        tbl_old[f"entity:{table_name}_id"] = list(map(lambda f: f['name'], ent_old))

        if store_membership:
            membership = list(map(lambda g: set(map(lambda h: h['entityName'], g['items'])), tbl_old['samples']))
            del tbl_old['samples']

        c = list(tbl_old.columns)
        c.remove(f"entity:{table_name}_id")
        c = [f"entity:{table_name}_id"] + c
        tbl_old = tbl_old[c]
        tbl_old = tbl_old.astype(str)

    return tbl_old, membership


def main():
    parser = argparse.ArgumentParser(description='Copy data to delivery workspaces', prog='deliver_data')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-i', '--copy-inputs', action='store_true', help="Copy input files too")
    args = parser.parse_args()

    tbl_old, _ = load_table(args.namespace, args.workspace, 'sample')
    ss_old, membership = load_table(args.namespace, args.workspace, 'sample_set', store_membership=True)

    tbl_filtered = tbl_old[~tbl_old.workspace.isin(['nan', ''])]

    print(f"Accessing Terra as '{fapi.whoami()}'")

    workspace_list = fapi.list_workspaces("workspace.name,workspace.namespace").json()
    workspaces = set()
    for w in workspace_list:
        workspaces.add(w['workspace']['name'])

    for index, row in tbl_filtered.iterrows():
        if row['workspace'] not in workspaces:
            a = fapi.create_workspace(args.namespace, row['workspace'])
            b = fapi.update_workspace_acl(args.namespace, row['workspace'], [
                {"email": "kiran@broadinstitute.org", "accessLevel": "OWNER"},
                {"email": "222581509023-compute@developer.gserviceaccount.com", "accessLevel": "OWNER"},
                #{"email": "shuang@broadinstitute.org", "accessLevel": "OWNER"},
                #{"email": "lholmes@broadinstitute.org", "accessLevel": "OWNER"},
            ])

            print(f"[workspace  : {a.status_code}] Created workspace '{row['workspace']}'")

        q = fapi.get_workspace(args.namespace, row['workspace']).json()

        if 'Garimella' in row['workspace']:
            newrow = row.replace('gs://broad-gp-pacbio-outgoing/', f"gs://{q['workspace']['bucketName']}/", regex=True)
            newrow.replace('gs://broad-gp-pacbio/', f"gs://{q['workspace']['bucketName']}/input/pacbio/", inplace=True, regex=True)
            newrow.replace('gs://broad-gp-oxfordnano-outgoing/', f"gs://{q['workspace']['bucketName']}/", inplace=True, regex=True)
            newrow.replace('gs://broad-gp-oxfordnano/', f"gs://{q['workspace']['bucketName']}/input/oxfordnano/", inplace=True, regex=True)

            a = fapi.copy_entities(args.namespace,
                                   args.workspace,
                                   args.namespace,
                                   newrow['workspace'],
                                   'sample',
                                   [newrow['entity:sample_id']])

            nr = newrow.to_dict()

            for k, v in row.to_dict().items():
                if 'gs://' in v:
                    if 'gs://broad-gp-pacbio/' in v or 'gs://broad-gp-oxfordnano/' in v:
                        if args.copy_inputs:
                            copy_file(v, nr[k])
                    else:
                        copy_file(v, nr[k])

            print(f"[sample     : {a.status_code}] Added '{row['entity:sample_id']}' to workspace '{row['workspace']}'")

            for i in range(len(membership)):
                if row['entity:sample_id'] in membership[i]:
                    a = fapi.copy_entities(args.namespace,
                                           args.workspace,
                                           args.namespace,
                                           newrow['workspace'],
                                           'sample_set',
                                           [ss_old['entity:sample_set_id'][i]],
                                           link_existing_entities=True)

                    print(f"[sample_set : {a.status_code}] Added '{ss_old['entity:sample_set_id'][i]}' to workspace '{row['workspace']}'")


if __name__ == "__main__":
    main()
