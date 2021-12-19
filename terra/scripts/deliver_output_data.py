#!/usr/bin/env python

import argparse
import subprocess
import re
import os

import pandas as pd
import firecloud.api as fapi
from google.cloud import storage

from multiprocessing.pool import Pool, ThreadPool
from functools import partial
from tqdm import tqdm


pd.set_option('max_columns', 200)
pd.set_option('max_rows', 200)
pd.set_option("max_colwidth", None)


def main():
    parser = argparse.ArgumentParser(description='Copy data to delivery workspaces', prog='deliver_output_data')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-u', '--user', type=str, nargs="+", help="Usernames to add as owners to the workspace")
    parser.add_argument('-t', '--threads', type=int, default=os.cpu_count() + 2, help="Number of threads to use to copy files")
    parser.add_argument('-r', '--run', action='store_true', help="Actually run the copying commands")
    args = parser.parse_args()

    print(f"Accessing Terra as '{fapi.whoami()}'.")
    print(f"Using {args.threads} threads to copy files.")

    owners = [{"email": u, "accessLevel": "OWNER"} for u in args.user]

    tbl_old, _ = load_table(args.namespace, args.workspace, 'sample')
    ss_old, membership = load_table(args.namespace, args.workspace, 'sample_set', store_membership=True)

    tbl_new_hash = {}
    ss_new_hash = {}
    membership_new_hash = {}
    namespace_new_hash = {}
    bucket_new_hash = {}
    copy_lists = {}

    tbl_filtered = tbl_old[~tbl_old.workspace.isin(['nan', ''])]
    workspaces = get_workspaces()

    for rw in tbl_filtered['workspace'].unique():
        if rw not in workspaces:
            create_workspace(args, rw)

            workspaces[rw] = args.namespace
            b = fapi.update_workspace_acl(args.namespace, rw, owners)

        if rw not in tbl_new_hash:
            tbl_new_hash[rw] = pd.DataFrame(columns=tbl_filtered.columns)
            ss_new_hash[rw] = pd.DataFrame(columns=ss_old.columns)
            membership_new_hash[rw] = []

        ns = workspaces[rw]

        tbl_rw = tbl_filtered[tbl_filtered['workspace'] == rw]
        qa = fapi.get_workspace(ns, rw)
        if qa.status_code == 200:
            q = qa.json()

            bucket_name = f"gs://{q['workspace']['bucketName']}"

            for index, row in tbl_rw.iterrows():
                newrow = row
                newrow = newrow.replace('gs://broad-gp-pacbio-outgoing/', bucket_name + "/", regex=True)
                newrow = newrow.replace('gs://broad-gp-oxfordnano-outgoing/', bucket_name + "/", regex=True)

                tbl_new_hash[rw] = tbl_new_hash[rw].append(newrow)
                bucket_new_hash[bucket_name] = rw

                for k, v in row.to_dict().items():
                    if 'gs://' in v and ('gs://broad-gp-pacbio/' not in v and 'gs://broad-gp-oxfordnano/' not in v):
                        if bucket_name not in copy_lists:
                            copy_lists[bucket_name] = {}

                        if '.' in v:
                            if len(re.sub("gs://", "", v).split("/")) <= 4:
                                copy_lists[bucket_name][v] = newrow[k]
                            else:
                                copy_lists[bucket_name][os.path.dirname(v) + "/"] = os.path.dirname(newrow[k]) + "/"

            for ss_index, ss_row in ss_old.iterrows():
                for sample_id in tbl_new_hash[rw]['entity:sample_id']:
                    if sample_id in membership[ss_index]:
                        new_ss_row = ss_row
                        new_ss_row = new_ss_row.replace('gs://broad-gp-pacbio-outgoing/', bucket_name + "/", regex=True)
                        new_ss_row = new_ss_row.replace('gs://broad-gp-oxfordnano-outgoing/', bucket_name + "/", regex=True)

                        for ks, vs in ss_row.to_dict().items():
                            if 'gs://' in vs and ('gs://broad-gp-pacbio/' not in vs and 'gs://broad-gp-oxfordnano/' not in vs):
                                if bucket_name not in copy_lists:
                                    copy_lists[bucket_name] = {}

                                if '.' in vs:
                                    if len(re.sub("gs://", "", vs).split("/")) <= 4:
                                        copy_lists[bucket_name][vs] = new_ss_row[ks]
                                    else:
                                        copy_lists[bucket_name][os.path.dirname(vs) + "/"] = os.path.dirname(new_ss_row[ks] + "/")

                        ss_new_hash[rw] = ss_new_hash[rw].append(ss_row)
                        membership_new_hash[rw].append(membership[ss_index])
                        namespace_new_hash[rw] = ns

    for workspace in membership_new_hash:
        if len(membership_new_hash[workspace]) > 0:
            oms = pd \
                .DataFrame({'entity:sample_set_id': list(ss_new_hash[workspace]['entity:sample_set_id']), 'sample': membership_new_hash[workspace]}) \
                .explode('sample', ignore_index=True)
            oms.columns = ['membership:sample_set_id', 'sample']

            if args.run:
                for ssname in list(oms['membership:sample_set_id']):
                    a = fapi.delete_sample_set(namespace_new_hash[workspace], workspace, ssname)

                upload_tables(namespace_new_hash[workspace], workspace, tbl_new_hash[workspace], ss_new_hash[workspace], oms)

    num_files = 0
    for bucket in copy_lists:
        for s in copy_lists[bucket]:
            num_files += 1

    print("Examining files...")
    num_files_copied = {}
    with tqdm(total=num_files) as pbar:
        pool = Pool(args.threads)
        for bucket in copy_lists:
            for s in copy_lists[bucket]:
                r = pool.apply_async(copy_files, args = (s, copy_lists[bucket][s], args.run, ))
                pbar.update()

                if bucket_new_hash[bucket] not in num_files_copied:
                    num_files_copied[bucket_new_hash[bucket]] = 0

                num_files_copied[bucket_new_hash[bucket]] += r.get()

        pool.close()
        pool.join()

    tot_files_copied = 0
    for workspace in num_files_copied:
        print(f"Workspace '{workspace}': copied {num_files_copied[workspace]} files. {'[dry-run]' if not args.run else ''}")
        tot_files_copied += num_files_copied[workspace]

    print(f"Total files copied: {tot_files_copied}. {'[dry-run]' if not args.run else ''}")

def create_workspace(args, rw, run):
    status_code = "000"
    if run:
        a = fapi.create_workspace(args.namespace, rw)
        status_code = a.status_code

    print(f"[workspace  : {status_code}] Created workspace '{rw}' {'[dry-run]' if run else ''}")


def copy_files(src, dst, run):
    if src[-1] == '/':
        if run:
            result = subprocess.run(["gsutil", "-m", "rsync", "-C", re.sub("/$", "", src), os.path.dirname(re.sub("/$", "", dst)) + "/"], capture_output=True, text=True)
            return len(list(filter(lambda x: 'Copying' in x, result.stderr.split("\n"))))
        else:
            result = subprocess.run(["gsutil", "-m", "rsync", "-C", "-n", re.sub("/$", "", src), os.path.dirname(re.sub("/$", "", dst)) + "/"], capture_output=True, text=True)
            return len(list(filter(lambda x: 'Would copy' in x, result.stderr.split("\n"))))
    else:
        if should_copy(src, dst):
            if run:
                result = subprocess.run(["gsutil", "cp", src, dst], capture_output=True, text=True)
            else:
                result = subprocess.run(["echo", "gsutil", "cp", src, dst, "[dry-run]"], capture_output=True, text=True)

            return 1

    return 0

def should_copy(src, dst):
    storage_client = storage.Client()

    bs, ns = re.sub("^gs://", "", src).split("/", maxsplit=1)
    bd, nd = re.sub("^gs://", "", dst).split("/", maxsplit=1)

    fs = storage.Blob(bucket=storage_client.bucket(bs), name=ns)
    if fs.exists():
        fs.reload()

    fd = storage.Blob(bucket=storage_client.bucket(bd), name=nd)
    if fd.exists():
        fd.reload()

    copy = not fd.exists() or fs.md5_hash != fd.md5_hash


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


def get_workspaces():
    workspace_list = fapi.list_workspaces("workspace.name,workspace.namespace").json()
    workspaces = {}
    for w in workspace_list:
        workspaces[w['workspace']['name']] = w['workspace']['namespace']
    return workspaces


if __name__ == "__main__":
    main()
