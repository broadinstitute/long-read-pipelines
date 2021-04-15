#!/usr/bin/python3

import os
import re
import hashlib

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


def upload_sample_set(namespace, workspace, tbl):
    # delete old sample set
    ss_old = fapi.get_entities(namespace, workspace, f'sample_set').json()
    sample_sets = list(map(lambda e: e['name'], ss_old))
    f = [fapi.delete_sample_set(namespace, workspace, sample_set_index) for sample_set_index in sample_sets]

    # upload new sample set
    ss = tbl.filter(['participant'], axis=1).drop_duplicates()
    ss.columns = [f'entity:sample_set_id']
    
    b = fapi.upload_entities(namespace, workspace, entity_data=ss.to_csv(index=False, sep="\t"), model='flexible')
    if b.status_code == 200:
        print(f'Uploaded {len(ss)} sample sets successfully.')
    else:
        print(b.json())
    
    # upload membership set
    ms = tbl.filter(['participant', 'entity:sample_id'], axis=1).drop_duplicates()
    ms.columns = [f'membership:sample_set_id', f'sample']
    
    c = fapi.upload_entities(namespace, workspace, entity_data=ms.to_csv(index=False, sep="\t"), model='flexible')
    if c.status_code == 200:
        print(f'Uploaded {len(ms)} sample set members successfully.')
    else:
        print(c.json())


namespace = 'broad-firecloud-dsde-methods'
workspace = 'dsp-oxfordnano'

gcs_buckets_ont = ['gs://broad-gp-oxfordnano', 'gs://broad-dsde-methods-long-reads-deepseq']

ent_old = fapi.get_entities(namespace, workspace, 'sample').json()

if len(ent_old) > 0:
    tbl_old = pd.DataFrame(list(map(lambda e: e['attributes'], ent_old)))
    tbl_old["entity:sample_id"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_old["final_summary_file"]))

ts = load_summaries(gcs_buckets_ont)

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

if len(ent_old) > 0:
    for sample_id in tbl_old.merge(tbl_new, how='outer', indicator=True).loc[lambda x : x['_merge']=='left_only']['entity:sample_id'].tolist():
        print(f'Entry for sample {sample_id} has been modified.  Keeping changes.')
        tbl_new = tbl_new[tbl_new['entity:sample_id'] != sample_id]

merged_tbl = pd.merge(tbl_old, tbl_new, how='outer') if len(ent_old) > 0 else tbl_new
merged_tbl = merged_tbl[['entity:sample_id'] + tbl_header]
merged_tbl = merged_tbl.drop_duplicates(subset=['flow_cell_id', 'instrument', 'sample_name', 'acquisition_stopped', 'processing_stopped', 'fast5_files_in_fallback', 'fast5_files_in_final_dest', 'fastq_files_in_fallback', 'fastq_files_in_final_dest'])
merged_tbl = merged_tbl[merged_tbl.fast5_files_in_final_dest != "0"]

a = fapi.upload_entities(namespace, workspace, entity_data=merged_tbl.to_csv(index=False, sep="\t"), model='flexible')

if a.status_code == 200:
    print(f'Uploaded {len(merged_tbl)} rows successfully.')
else:
    print(a.json())

upload_sample_set(namespace, workspace, merged_tbl)
