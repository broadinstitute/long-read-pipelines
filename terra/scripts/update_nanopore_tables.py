#!/usr/bin/env python

import os
import re
import math
import hashlib
import argparse

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


def load_summaries(gcs_buckets, project):
    storage_client = storage.Client(project=project)

    ts = []
    for gcs_bucket in gcs_buckets:
        blobs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket))

        for blob in blobs:
            if 'final_summary' in blob.name:
                doc = blob.download_as_string()
                t = {}

                for line in doc.decode("utf-8").split("\n"):
                    if '=' in line:
                        k, v = line.split('=')

                        t[k] = v

                t['Files'] = {
                    'fast5_pass_dir': gcs_bucket + "/" + os.path.dirname(blob.name) + "/fast5_pass",
                    'fastq_pass_dir': gcs_bucket + "/" + os.path.dirname(blob.name) + "/fastq_pass",
                    'final_summary.txt': gcs_bucket + "/" + blob.name,
                    'sequencing_summary.txt': 'missing'
                }

                bs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket),
                                               prefix=os.path.dirname(blob.name) + "/" + t['sequencing_summary_file'])
                for b in bs:
                    t['Files']['sequencing_summary.txt'] = gcs_bucket + "/" + b.name

                if 'sequencing_summary.txt' not in t['Files']:
                    pp = pprint.PrettyPrinter(indent=4)
                    pp.pprint(t)

                ts.append(t)

    return ts


def load_new_sample_table(buckets, project):
    ts = load_summaries(buckets, project)

    tbl_header = ["final_summary_file", "sequencing_summary_file", "fast5_pass_dir", "fastq_pass_dir", "protocol_group_id", "instrument", "position", "flow_cell_id", "sample_name", "basecalling_enabled", "started", "acquisition_stopped", "processing_stopped", "fast5_files_in_fallback", "fast5_files_in_final_dest", "fastq_files_in_fallback", "fastq_files_in_final_dest"]
    tbl_rows = []

    for e in ts:
        tbl_rows.append([
            e["Files"]["final_summary.txt"],
            e["Files"]["sequencing_summary.txt"],
            e["Files"]["fast5_pass_dir"],
            e["Files"]["fastq_pass_dir"],
            e["protocol_group_id"],
            e["instrument"],
            e["position"],
            e["flow_cell_id"],
            e["sample_id"],
            e["basecalling_enabled"],
            e["started"],
            e["acquisition_stopped"],
            e["processing_stopped"],
            e["fast5_files_in_fallback"],
            e["fast5_files_in_final_dest"],
            e["fastq_files_in_fallback"],
            e["fastq_files_in_final_dest"]
        ])

    tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)
    tbl_new["entity:sample_id"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_new["final_summary_file"]))
    tbl_new = tbl_new.astype(str)

    return tbl_new


def merge_tables(tbl_old, tbl_new):
    if tbl_old is not None:
        outer_tbl = pd.merge(tbl_old, tbl_new, how='outer', sort=True, indicator=True)
    else:
        outer_tbl = tbl_new

    hs = []
    for l in list(outer_tbl['entity:sample_id'].unique()):
        g = outer_tbl.loc[outer_tbl['entity:sample_id'] == l].sort_values('_merge')

        if len(g) == 1:
            hs.append(g.iloc[0].to_dict())
        else:
            h = {}
            for col_name in list(outer_tbl.columns):
                q = g[col_name]
                v = q.where((q != 'None') & (q != 'nan')).dropna()
                h[col_name] = v.iloc[0] if len(v) > 0 else ''

            hs.append(h)

    joined_tbl = pd.DataFrame(hs)

    if '_merge' in joined_tbl:
        del joined_tbl['_merge']
    c = list(joined_tbl.columns)
    c.remove("entity:sample_id")
    c = ["entity:sample_id"] + c
    joined_tbl = joined_tbl[c]

    joined_tbl['sample_name'] = joined_tbl['sample_name'].str.replace(r'\s+', ' ', regex=True).astype('str')

    return joined_tbl


def update_sample_table(namespace, workspace, buckets, project):
    tbl_old, _ = load_table(namespace, workspace, 'sample')
    tbl_new = load_new_sample_table(buckets, project)
    joined_tbl = merge_tables(tbl_old, tbl_new)
    joined_tbl = joined_tbl.replace('^nan$', '', regex=True)

    return joined_tbl


def update_sample_set_table(namespace, workspace, joined_tbl):
    ss_old, membership = load_table(namespace, workspace, 'sample_set', store_membership=True)

    # create old membership set
    oms = pd \
        .DataFrame({'entity:sample_set_id': list(ss_old['entity:sample_set_id']), 'sample': membership}) \
        .explode('sample', ignore_index=True)
    oms.columns = ['membership:sample_set_id', 'sample']

    # create sample set
    ss = joined_tbl.filter(['sample_name'], axis=1).drop_duplicates()
    ss.columns = [f'entity:sample_set_id']

    if ss_old is not None:
        ss = pd.merge(ss_old, ss, how='outer', sort=True)
    ss = ss.replace('nan', '', regex=True)

    # create new membership set
    ms = joined_tbl.filter(['sample_name', 'entity:sample_id'], axis=1).drop_duplicates()
    ms.columns = [f'membership:sample_set_id', f'sample']

    # create full membership set
    fms = pd.merge(ms, oms, how='outer', indicator=True)

    # create new/modified membership set
    nms = fms[fms['_merge'] != 'both']

    return ss, nms


def upload_table(namespace, workspace, table, label):
    # upload new samples
    a = fapi.upload_entities(namespace, workspace, entity_data=table.to_csv(index=False, sep="\t"), model='flexible')

    if a.status_code == 200:
        print(f'Uploaded {len(table)} {label} rows successfully.')
    else:
        print(a.json())


def upload_tables(namespace, workspace, s, ss, nms):
    for ssname in list(nms[nms['_merge'] == 'right_only']['membership:sample_set_id']):
        a = fapi.delete_sample_set(namespace, workspace, ssname)

        if a.status_code == 204:
            print(f'Removed out-of-date sample set {ssname} successfully.')
        else:
            print(a.json())

    lms = nms[nms['_merge'] == 'left_only'][['membership:sample_set_id', 'sample']]

    upload_table(namespace, workspace, s, 'sample')
    upload_table(namespace, workspace, ss, 'sample_set')
    upload_table(namespace, workspace, lms, 'sample_set membership')


def main():
    parser = argparse.ArgumentParser(description='Update Terra workspace sample table', prog='update_nanopore_tables')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-r', '--run', action='store_true', help="Turn off the default dry-run mode")
    parser.add_argument('buckets', metavar='B', type=str, nargs='+', help='GCS buckets to scan')
    args = parser.parse_args()

    s = update_sample_table(args.namespace, args.workspace, args.buckets, args.project)
    ss, nms = update_sample_set_table(args.namespace, args.workspace, s)

    if args.run:
        upload_tables(args.namespace, args.workspace, s, ss, nms)


if __name__ == "__main__":
    main()
