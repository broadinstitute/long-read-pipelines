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


def traverse_xml(key, xml):
    tables = []
    table = {}

    for k in xml:
        if 'xmlns' in k or 'xsi' in k:
            continue

        v = xml[k]

        k = re.sub('^@|^#|^pbds:|^pbbase:', '', k)

        l = []
        if isinstance(v, str) or isinstance(v, dict):
            l = [v]
        elif isinstance(v, list):
            l = v

        for va in l:
            if isinstance(va, str):
                table[k] = v
            if isinstance(va, dict):
                f = traverse_xml(k, va)
                tables.extend(f)

    if len(table) > 0:
        tables.append({key: table})

    return tables


def combine(tables):
    combined_tables = {}

    for table in tables:
        for k in table:
            if k not in combined_tables:
                combined_tables[k] = []

            combined_tables[k].append(table[k])

    return combined_tables


def load_xmls(gcs_buckets):
    storage_client = storage.Client(project='broad-dsp-lrma')
    schemas = OrderedDict()

    ts = []
    for gcs_bucket in gcs_buckets:
        blobs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket))

        for blob in blobs:
            if 'subreadset.xml' in blob.name:
                xml = blob.download_as_string()
                doc = xmltodict.parse(xml)

                t = combine(traverse_xml('root', doc))
                t['Files'] = {
                    'subreadset.xml': gcs_bucket + "/" + blob.name,
                    'subreads.bam': gcs_bucket + "/" + re.sub("et.xml", ".bam", blob.name),
                    'subreads.bam.pbi': gcs_bucket + "/" + re.sub("et.xml", ".bam.pbi", blob.name)
                }
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
workspace = 'dsp-pacbio'
gcs_buckets_pb = ['gs://broad-gp-pacbio']

ent_old = fapi.get_entities(namespace, workspace, 'sample').json()

if len(ent_old) > 0:
    tbl_old = pd.DataFrame(list(map(lambda e: e['attributes'], ent_old)))
    tbl_old["entity:sample_id"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_old["subreads_bam"]))

ts = load_xmls(gcs_buckets_pb)

tbl_header = ["subreads_bam", "subreads_pbi", "movie_name", "well_name", "created_at", "original_participant_name", "participant", "insert_size", "is_ccs", "isoseq"]
tbl_rows = []

for e in ts:
    tbl_rows.append([
        e['Files']['subreads.bam'],
        e['Files']['subreads.bam.pbi'],
        e['CollectionMetadata'][0]['Context'] if 'Context' in e['CollectionMetadata'][0] else "UnknownFlowcell",
        e['WellSample'][0]['WellName'] if 'WellName' in e['WellSample'][0] else "Z00",
        e['WellSample'][0]['CreatedAt'] if 'CreatedAt' in e['WellSample'][0] else "0001-01-01T00:00:00",
        re.sub("[# ]", "", e['WellSample'][0]['Name']) if 'Name' in e['WellSample'][0] else "UnknownSample",
        re.sub("[# ]", "", e['WellSample'][0]['Name']) if 'Name' in e['WellSample'][0] else "UnknownSample",        
        e['WellSample'][0]['InsertSize'] if 'InsertSize' in e['WellSample'][0] else "0",
        e['WellSample'][0]['IsCCS'] if 'IsCCS' in e['WellSample'][0] else "unknown",
        e['WellSample'][0]['IsoSeq'] if 'IsoSeq' in e['WellSample'][0] else "unknown"
    ])
    
tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)
tbl_new["entity:sample_id"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_new["subreads_bam"]))
tbl_new["entity:ccs_sample"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_new["subreads_bam"]))
tbl_new["entity:clr_sample"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_new["subreads_bam"]))
tbl_new["entity:isoseq_sample"] = list(map(lambda f: hashlib.md5(f.encode("utf-8")).hexdigest(), tbl_new["subreads_bam"]))

if len(ent_old) > 0:
    for sample_id in tbl_old.merge(tbl_new, how='outer', indicator=True).loc[lambda x : x['_merge']=='left_only']['entity:sample_id'].tolist():
        print(f'Entry for sample {sample_id} has been modified.  Keeping changes.')
        tbl_new = tbl_new[tbl_new['entity:sample_id'] != sample_id]

merged_tbl = pd.merge(tbl_old, tbl_new, how='outer') if len(ent_old) > 0 else tbl_new
merged_tbl = merged_tbl[['entity:sample_id'] + tbl_header]

a = fapi.upload_entities(namespace, workspace, entity_data=merged_tbl.to_csv(index=False, sep="\t"), model='flexible')

if a.status_code == 200:
    print(f'Uploaded {len(merged_tbl)} rows successfully.')
else:
    print(a.json())

upload_sample_set(namespace, workspace, merged_tbl)
