#!/usr/bin/env python

import argparse
import subprocess
import re
import os

import pandas as pd
import firecloud.api as fapi

from multiprocessing.pool import Pool, ThreadPool
from functools import partial
from tqdm import tqdm


pd.set_option('max_columns', 200)
pd.set_option('max_rows', 200)
pd.set_option("max_colwidth", None)


def copy_file(src, dst, run):
    if run:
        result = subprocess.run(["gsutil", "cp", "-n", src, dst], capture_output=True, text=True)
        return "Skipping" not in result.stdout
    else:
        result = subprocess.run(["echo", "gsutil", "cp", "-n", src, dst, "[dry-run]"], capture_output=True, text=True)

    return False


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
    parser = argparse.ArgumentParser(description='Copy input data to delivery workspaces', prog='deliver_input_data')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-t', '--threads', type=int, default=os.cpu_count() + 2, help="Number of threads to use to copy files")
    parser.add_argument('-r', '--run', action='store_true', help="Actually run the copying commands")
    args = parser.parse_args()

    tbl_old, _ = load_table(args.namespace, args.workspace, 'sample')
    ss_old, membership = load_table(args.namespace, args.workspace, 'sample_set', store_membership=True)

    tbl_filtered = tbl_old[~tbl_old.workspace.isin(['nan', ''])]

    print(f"Accessing Terra as '{fapi.whoami()}'.")
    print(f"Using {args.threads} threads to copy files.")

    tbl_new = pd.DataFrame(columns=tbl_filtered.columns)
    ss_new = pd.DataFrame(columns=ss_old.columns)
    membership_new = []

    copy_lists = {}

    for index, row in tbl_filtered.iterrows():
        qa = fapi.get_workspace(args.namespace, args.workspace)

        if qa.status_code == 200:
            q = qa.json()

            bucket_name = f"gs://{q['workspace']['bucketName']}"
            newrow = row.replace('gs://broad-gp-pacbio/', bucket_name + "/inputs/pacbio/", regex=True)
            newrow.replace('gs://broad-gp-oxfordnano/', bucket_name + "/inputs/oxfordnano/", inplace=True, regex=True)

            tbl_new = tbl_new.append(newrow)

            for k, v in row.to_dict().items():
                if 'gs://' in v:
                    if 'gs://broad-gp-pacbio/' in v or 'gs://broad-gp-oxfordnano/' in v:
                        if bucket_name not in copy_lists:
                            copy_lists[bucket_name] = {}
                        if '.' in v:
                            copy_lists[bucket_name][v] = newrow[k]

            for ss_index, ss_row in ss_old.iterrows():
                if row['entity:sample_id'] in membership[ss_index]:
                    ss_new = ss_new.append(ss_row)
                    membership_new.append(membership[ss_index])

    if len(membership_new) > 0:
        oms = pd \
            .DataFrame({'entity:sample_set_id': list(ss_new['entity:sample_set_id']), 'sample': membership_new}) \
            .explode('sample', ignore_index=True)
        oms.columns = ['membership:sample_set_id', 'sample']

        if args.run or True:
            for ssname in list(oms['membership:sample_set_id']):
                a = fapi.delete_sample_set(args.namespace, args.workspace, ssname)

            upload_tables(args.namespace, args.workspace, tbl_new, ss_new, oms)

    num_files = 0
    for bucket in copy_lists:
        for s in copy_lists[bucket]:
            num_files += 1

    print("Examining files...")
    num_files_copied = 0
    with tqdm(total=num_files) as pbar:
        pool = Pool(args.threads)
        for bucket in copy_lists:
            for s in copy_lists[bucket]:
                r = pool.apply_async(copy_file, args = (s, copy_lists[bucket][s], args.run, ))
                pbar.update()

                if r.get():
                    num_files_copied += 1

        pool.close()
        pool.join()

    print(f"Workspace '{args.workspace}': copied {num_files_copied} raw input files. {'[dry-run]' if not args.run else ''}")


if __name__ == "__main__":
    main()
