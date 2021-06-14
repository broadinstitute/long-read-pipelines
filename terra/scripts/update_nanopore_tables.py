#!/usr/bin/python3

import os
import re
import hashlib
import argparse

import pandas as pd
import firecloud.api as fapi

from google.cloud import bigquery
from google.cloud import storage
from google.api_core.exceptions import NotFound

from collections import OrderedDict

import xmltodict
import pprint


def load_summaries(gcs_buckets):
    storage_client = storage.Client(project='broad-dsp-lrma')
    schemas = OrderedDict()

    ts = []
    for gcs_bucket in gcs_buckets:
        blobs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket))

        for blob in blobs:
            if 'final_summary' in blob.name:
                doc = blob.download_as_string()
                t = {}
                
                for line in doc.decode("utf-8").split("\n"):
                    if '=' in line:
                        k,v = line.split('=')
                    
                        t[k] = v
                        
                t['Files'] = {
                    'final_summary.txt': gcs_bucket + "/" + blob.name,
                    'sequencing_summary.txt': 'missing'
                }

                bs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket), prefix=os.path.dirname(blob.name) + "/" + t['sequencing_summary_file'])
                for b in bs:
                    t['Files']['sequencing_summary.txt'] = gcs_bucket + "/" + b.name
                    
                if 'sequencing_summary.txt' not in t['Files']:
                    pp = pprint.PrettyPrinter(indent=4)
                    pp.pprint(t)

                ts.append(t)

    return ts


def upload_data(namespace, workspace, tbl):
    # delete old sample set
    ss_old = fapi.get_entities(namespace, workspace, f'sample_set').json()
    sample_sets = list(map(lambda e: e['name'], ss_old))
    f = [fapi.delete_sample_set(namespace, workspace, sample_set_index) for sample_set_index in sample_sets]
    
    # delete old samples
    s_old = fapi.get_entities(namespace, workspace, 'sample').json()
    samples = list(map(lambda e: e['name'], s_old))
    f = [fapi.delete_sample(namespace, workspace, sample_index) for sample_index in samples]

    # upload new samples
    a = fapi.upload_entities(namespace, workspace, entity_data=tbl.to_csv(index=False, sep="\t"), model='flexible')

    if a.status_code == 200:
        print(f'Uploaded {len(tbl)} rows successfully.')
    else:
        print(a.json())

    # upload new sample set
    ss = tbl.filter(['sample_name'], axis=1).drop_duplicates()
    ss.columns = [f'entity:sample_set_id']
    
    b = fapi.upload_entities(namespace, workspace, entity_data=ss.to_csv(index=False, sep="\t"), model='flexible')
    if b.status_code == 200:
        print(f'Uploaded {len(ss)} sample sets successfully.')
    else:
        print(b.json())
    
    # upload membership set
    ms = tbl.filter(['sample_name', 'entity:sample_id'], axis=1).drop_duplicates()
    ms.columns = [f'membership:sample_set_id', f'sample']
    
    c = fapi.upload_entities(namespace, workspace, entity_data=ms.to_csv(index=False, sep="\t"), model='flexible')
    if c.status_code == 200:
        print(f'Uploaded {len(ms)} sample set members successfully.')
    else:
        print(c.json())


def main():
    parser = argparse.ArgumentParser(description='Update Terra workspace sample table', prog='update_nanopore_tables')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('buckets', metavar='B', type=str, nargs='+', help='GCS buckets to scan')
    args = parser.parse_args()

    ent_old = fapi.get_entities(args.namespace, args.workspace, 'sample').json()
    tbl_old = None

    if len(ent_old) > 0:
        tbl_old = pd.DataFrame(list(map(lambda e: e['attributes'], ent_old)))
        tbl_old["entity:sample_id"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_old["final_summary_file"]))

    ts = load_summaries(args.buckets)

    tbl_header = ["final_summary_file", "sequencing_summary_file", "protocol_group_id", "instrument", "position", "flow_cell_id", "sample_name", "basecalling_enabled", "started", "acquisition_stopped", "processing_stopped", "fast5_files_in_fallback", "fast5_files_in_final_dest", "fastq_files_in_fallback", "fastq_files_in_final_dest"]
    tbl_rows = []

    for e in ts:
        tbl_rows.append([
            e["Files"]["final_summary.txt"],
            e["Files"]["sequencing_summary.txt"],
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

    if tbl_old is not None:
        outer_tbl = pd.merge(tbl_old, tbl_new, how='outer', sort=True, indicator=True)
    else:
        outer_tbl = tbl_new

    hs = []
    for l in list(outer_tbl['entity:sample_id'].unique()):
        g = outer_tbl.loc[outer_tbl['entity:sample_id'] == l]

        if len(g) == 1:
            hs.append(g.iloc[0].to_dict())
        else:
            h = {}
            for col_name in list(outer_tbl.columns):
                side = "left_only" if col_name in list(tbl_old.columns) else "right_only"
                q = list(g.loc[g['_merge'] == side][col_name])
                if len(q) > 0:
                    h[col_name] = q[0]

            hs.append(h)

    joined_tbl = pd.DataFrame(hs)
    if '_merge' in joined_tbl:
        del joined_tbl['_merge']

    c = list(joined_tbl.columns)
    c.remove("entity:sample_id")
    c = ["entity:sample_id"] + c
    joined_tbl = joined_tbl[c]

    upload_data(args.namespace, args.workspace, joined_tbl)


if __name__ == "__main__":
    main()
