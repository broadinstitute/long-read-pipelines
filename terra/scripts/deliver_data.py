#!/usr/bin/env python

import argparse
import subprocess
import re
import os

import pandas as pd
import firecloud.api as fapi

from multiprocessing.pool import Pool, ThreadPool
from functools import partial


pd.set_option('max_columns', 200)
pd.set_option('max_rows', 200)
pd.set_option("max_colwidth", None)


def copy_file(src, dst, run):
    if run:
        subprocess.run(["gsutil", "cp", "-n", src, dst])
    else:
        subprocess.run(["echo", "gsutil", "cp", "-n", src, dst, "[dry-run]"])


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


def upload_table(namespace, workspace, table, label):
    # upload new samples
    a = fapi.upload_entities(namespace, workspace, entity_data=table.to_csv(index=False, sep="\t"), model='flexible')

    if a.status_code == 200:
        print(f'Uploaded {len(table)} {label} rows successfully.')
    else:
        print(a.json())


def upload_tables(namespace, workspace, s, ss, nms):
    upload_table(namespace, workspace, s, 'sample')
    upload_table(namespace, workspace, ss, 'sample_set')
    upload_table(namespace, workspace, nms, 'sample_set membership')


def main():
    parser = argparse.ArgumentParser(description='Copy data to delivery workspaces', prog='deliver_data')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-u', '--user', type=str, nargs="+", help="Usernames to add as owners to the workspace")
    parser.add_argument('-i', '--copy-inputs', action='store_true', help="Copy input files too")
    parser.add_argument('-r', '--run', action='store_true', help="Actually run the copying commands")
    args = parser.parse_args()

    tbl_old, _ = load_table(args.namespace, args.workspace, 'sample')
    ss_old, membership = load_table(args.namespace, args.workspace, 'sample_set', store_membership=True)

    tbl_filtered = tbl_old[~tbl_old.workspace.isin(['nan', ''])]

    print(f"Accessing Terra as '{fapi.whoami()}'")

    workspace_list = fapi.list_workspaces("workspace.name,workspace.namespace").json()
    workspaces = {}
    for w in workspace_list:
        workspaces[w['workspace']['name']] = w['workspace']['namespace']

    tbl_new_hash = {}
    ss_new_hash = {}
    membership_new_hash = {}
    namespace_new_hash = {}

    copy_lists = {}

    for index, row in tbl_filtered.iterrows():
        rw = row['workspace']

        if rw not in workspaces:
            owners = [{"email": u, "accessLevel": "OWNER"} for u in args.user]

            if args.run:
                a = fapi.create_workspace(args.namespace, rw)
                b = fapi.update_workspace_acl(args.namespace, rw, owners)

                # b = fapi.update_workspace_acl(args.namespace, rw, [
                #     {"email": "kiran@broadinstitute.org", "accessLevel": "OWNER"},
                #     {"email": "222581509023-compute@developer.gserviceaccount.com", "accessLevel": "OWNER"},
                #     {"email": "shuang@broadinstitute.org", "accessLevel": "OWNER"},
                #     {"email": "lholmes@broadinstitute.org", "accessLevel": "OWNER"},
                # ])

                print(f"[workspace  : {a.status_code}] Created workspace '{rw}'")
            else:
                print(f"[workspace  : 000] Created workspace '{rw}' [dry-run]")

            workspaces[rw] = args.namespace

        ns = workspaces[rw]

        if rw not in tbl_new_hash:
            tbl_new_hash[rw] = pd.DataFrame(columns=tbl_filtered.columns)
            ss_new_hash[rw] = pd.DataFrame(columns=ss_old.columns)
            membership_new_hash[rw] = []

        qa = fapi.get_workspace(ns, rw)

        if qa.status_code == 200:
            q = qa.json()

            bucket_name = f"gs://{q['workspace']['bucketName']}"
            newrow = row.replace('gs://broad-gp-pacbio-outgoing/', bucket_name + "/", regex=True)
            newrow.replace('gs://broad-gp-pacbio/', bucket_name + "/inputs/pacbio/", inplace=True, regex=True)
            newrow.replace('gs://broad-gp-oxfordnano-outgoing/', bucket_name + "/", inplace=True, regex=True)
            newrow.replace('gs://broad-gp-oxfordnano/', bucket_name + "/inputs/oxfordnano/", inplace=True, regex=True)

            tbl_new_hash[rw] = tbl_new_hash[rw].append(newrow)

            for k, v in row.to_dict().items():
                if 'gs://' in v:
                    if 'gs://broad-gp-pacbio/' in v or 'gs://broad-gp-oxfordnano/' in v:
                        if args.copy_inputs:
                            if bucket_name not in copy_lists:
                                copy_lists[bucket_name] = {}
                            if '.' in v:
                                copy_lists[bucket_name][v] = newrow[k]
                    else:
                        if bucket_name not in copy_lists:
                            copy_lists[bucket_name] = {}
                        if '.' in v:
                            copy_lists[bucket_name][v] = newrow[k]

            for ss_index, ss_row in ss_old.iterrows():
                if row['entity:sample_id'] in membership[ss_index]:
                    ss_new_hash[rw] = ss_new_hash[rw].append(ss_row)
                    membership_new_hash[rw].append(membership[ss_index])
                    namespace_new_hash[rw] = ns

    for workspace in membership_new_hash:
        print(f'Processing workspace {workspace}...')

        if len(membership_new_hash[workspace]) > 0:
            oms = pd \
                .DataFrame({'entity:sample_set_id': list(ss_new_hash[workspace]['entity:sample_set_id']), 'sample': membership_new_hash[workspace]}) \
                .explode('sample', ignore_index=True)
            oms.columns = ['membership:sample_set_id', 'sample']

            if args.run:
                for ssname in list(oms['membership:sample_set_id']):
                    a = fapi.delete_sample_set(namespace_new_hash[workspace], workspace, ssname)

                upload_tables(namespace_new_hash[workspace], workspace, tbl_new_hash[workspace], ss_new_hash[workspace], oms)
            else:
                print(f'\tUploaded {len(tbl_new_hash[workspace])} rows successfully. [dry-run]')


    num_files = 0
    num_threads = os.cpu_count() + 2
    pool = Pool(num_threads)
    for bucket in copy_lists:
        for s in copy_lists[bucket]:
            pool.apply_async(copy_file, args = (s, copy_lists[bucket][s], args.run, ))
            num_files += 1

    pool.close()
    pool.join()

    print(f"Copied {num_files} files using {num_threads} threads.")


if __name__ == "__main__":
    main()
