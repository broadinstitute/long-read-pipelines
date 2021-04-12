#!/usr/bin/env python

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

        k = re.sub('^@|^#|^pbds:|^pbbase:|^pbmeta:|^pbsample:', '', k)

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
    storage_client = storage.Client()
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
                    'input_dir': os.path.dirname(gcs_bucket + "/" + blob.name),
                    
                    'subreadset.xml': gcs_bucket + "/" + blob.name,
                    'subreads.bam': gcs_bucket + "/" + re.sub("et.xml", ".bam", blob.name),
                    'subreads.bam.pbi': gcs_bucket + "/" + re.sub("et.xml", ".bam.pbi", blob.name),
                    
                    'consensusreadset.xml': "",
                    'ccs_reports.txt': "",
                    'reads.bam': "",
                    'reads.bam.pbi': ""
                }
                ts.append(t)
            elif 'consensusreadset.xml' in blob.name:
                xml = blob.download_as_string()
                doc = xmltodict.parse(xml)

                t = combine(traverse_xml('root', doc))
                t['Files'] = {
                    'input_dir': os.path.dirname(gcs_bucket + "/" + blob.name),
                    
                    'subreadset.xml': "",
                    'subreads.bam': "",
                    'subreads.bam.pbi': "",

                    'consensusreadset.xml': gcs_bucket + "/" + blob.name,
                    'ccs_reports.txt': gcs_bucket + "/" + re.sub(".consensusreadset.xml", ".ccs_reports.txt", blob.name),
                    'reads.bam': gcs_bucket + "/" + re.sub(".consensusreadset.xml", ".reads.bam", blob.name),
                    'reads.bam.pbi': gcs_bucket + "/" + re.sub(".consensusreadset.xml", ".reads.bam.pbi", blob.name)
                }
                ts.append(t)

    return ts


def upload_samples(namespace, workspace, tbl):
    a = fapi.upload_entities(namespace, workspace, entity_data=tbl.to_csv(index=False, sep="\t"), model='flexible')

    if a.status_code == 200:
        print(f'Uploaded {len(tbl)} rows successfully.')
    else:
        print(a.json())


def upload_sample_set(namespace, workspace, tbl):
    # delete old sample set
    ss_old = fapi.get_entities(namespace, workspace, f'sample_set').json()
    sample_sets = list(map(lambda e: e['name'], ss_old))
    f = [fapi.delete_sample_set(namespace, workspace, sample_set_index) for sample_set_index in sample_sets]

    # upload new sample set
    ss = tbl.filter(['bio_sample'], axis=1).drop_duplicates()
    ss.columns = [f'entity:sample_set_id']
    
    b = fapi.upload_entities(namespace, workspace, entity_data=ss.to_csv(index=False, sep="\t"), model='flexible')
    if b.status_code == 200:
        print(f'Uploaded {len(ss)} sample sets successfully.')
    else:
        print(b.json())
    
    # upload membership set
    ms = tbl.filter(['bio_sample', 'entity:sample_id'], axis=1).drop_duplicates()
    ms.columns = [f'membership:sample_set_id', f'sample']
    
    c = fapi.upload_entities(namespace, workspace, entity_data=ms.to_csv(index=False, sep="\t"), model='flexible')
    if c.status_code == 200:
        print(f'Uploaded {len(ms)} sample set members successfully.')
    else:
        print(c.json())


def load_ccs_report(ccs_report_path):
    d = {
        'ZMWs input': "",
        'ZMWs pass filters': "",
        'ZMWs fail filters': "",
        'ZMWs shortcut filters': "",
        'ZMWs with tandem repeats': "",
        'Below SNR threshold': "",
        'Median length filter': "",
        'Lacking full passes': "",
        'Heteroduplex insertions': "",
        'Coverage drops': "",
        'Insufficient draft cov': "",
        'Draft too different': "",
        'Draft generation error': "",
        'Draft above --max-length': "",
        'Draft below --min-length': "",
        'Reads failed polishing': "",
        'Empty coverage windows': "",
        'CCS did not converge': "",
        'CCS below minimum RQ': "",
        'Unknown error': ""
    }
    
    if ccs_report_path != "":
        storage_client = storage.Client()

        ccs_report = re.sub("^gs://", "", e['Files']['ccs_reports.txt']).split("/")
        blobs = storage_client.list_blobs(ccs_report[0], prefix="/".join(ccs_report[1:]))

        for blob in blobs:
            blob.download_to_filename("ccs_report.txt")

            file = open("ccs_report.txt", "r")

            d = {}
            for line in file:
                if len(line) > 1 and 'Exclusive counts for ZMWs' not in line:
                    a = line.rstrip().split(":")

                    k = a[0].rstrip()
                    v = float(re.sub(" ", "", re.sub(" \(.*$", "", a[1])))

                    if k not in d:
                        d[k] = 0.0;

                    d[k] = d[k] + v

            break
            
    return d


namespace = "broad-firecloud-dsde-methods" 
workspace = "dsp-pacbio"
gcs_buckets_pb = ['gs://broad-gp-pacbio']

ent_old = fapi.get_entities(namespace, workspace, 'sample').json()

if len(ent_old) > 0:
    tbl_old = pd.DataFrame(list(map(lambda e: e['attributes'], ent_old)))
    tbl_old["entity:sample_id"] = list(map(lambda f: f['name'], ent_old))

ts = load_xmls(gcs_buckets_pb)

tbl_header = ["entity:sample_id", "instrument", "movie_name", "well_name", "created_at", "bio_sample", "well_sample", "insert_size", "is_ccs", "is_isoseq", "num_records", "total_length", "zmws_input", "zmws_pass", "zmws_fail", "zmws_shortcut_filters", "gcs_input_dir", "subreads_bam", "subreads_pbi", "ccs_bam", "ccs_pbi"]
tbl_rows = []

for e in ts:
    r = load_ccs_report(e['Files']['ccs_reports.txt'])
    
    tbl_rows.append([
        e['CollectionMetadata'][0]['UniqueId'] if 'Context' in e['CollectionMetadata'][0] else "",

        e['CollectionMetadata'][0]['InstrumentName'] if 'Context' in e['CollectionMetadata'][0] else "UnknownInstrument",
        e['CollectionMetadata'][0]['Context'] if 'Context' in e['CollectionMetadata'][0] else "UnknownFlowcell",
        e['WellSample'][0]['WellName'] if 'WellName' in e['WellSample'][0] else "Z00",
        e['WellSample'][0]['CreatedAt'] if 'CreatedAt' in e['WellSample'][0] else "0001-01-01T00:00:00",
        re.sub("[# ]", "", e['BioSample'][0]['Name']) if 'BioSample' in e else "UnknownBioSample",
        re.sub("[# ]", "", e['WellSample'][0]['Name']) if 'Name' in e['WellSample'][0] else "UnknownWellSample",
        e['WellSample'][0]['InsertSize'] if 'InsertSize' in e['WellSample'][0] else "0",
        e['WellSample'][0]['IsCCS'] if 'IsCCS' in e['WellSample'][0] else "unknown",
        e['WellSample'][0]['IsoSeq'] if 'IsoSeq' in e['WellSample'][0] else "unknown",
        
        e['DataSetMetadata'][0]['NumRecords'],
        e['DataSetMetadata'][0]['TotalLength'],
        
        r['ZMWs input'],
        r['ZMWs pass filters'],
        r['ZMWs fail filters'],
        r['ZMWs shortcut filters'],

        e['Files']['input_dir'],
        e['Files']['subreads.bam'],
        e['Files']['subreads.bam.pbi'],
        e['Files']['reads.bam'],
        e['Files']['reads.bam.pbi'],
    ])
    
tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)

if len(ent_old) > 0:
    for sample_id in tbl_old.merge(tbl_new, how='outer', indicator=True).loc[lambda x : x['_merge']=='left_only']['entity:sample_id'].tolist():
        print(f'Entry for sample {sample_id} has been modified.  Keeping changes.')
        tbl_new = tbl_new[tbl_new['entity:sample_id'] != sample_id]
        
merged_tbl = pd.merge(tbl_old, tbl_new, how='outer') if len(ent_old) > 0 else tbl_new
merged_tbl = merged_tbl.drop_duplicates(subset=['entity:sample_id'])

l = list(merged_tbl.columns)
l.remove("entity:sample_id")
l = ["entity:sample_id"] + l
merged_tbl = merged_tbl[l]

upload_samples(namespace, workspace, merged_tbl)
upload_sample_set(namespace, workspace, merged_tbl)
