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

# from google.cloud import bigquery
from google.cloud import storage
from google.api_core.exceptions import NotFound

from collections import OrderedDict

from tqdm import tqdm
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


def load_new_sample_table(buckets, project):
    ts = load_xmls(buckets, project)
    tbl_header = ["entity:sample_id", "flowcell_id", "instrument", "movie_name", "well_name", "created_at", "bio_sample",
                  "well_sample", "insert_size", "is_ccs", "is_isoseq", "is_corrected", "description", "application",
                  "experiment_type", "num_records", "total_length", "ccs_report", "ccs_zmws_input",
                  "ccs_zmws_pass_filters", "ccs_zmws_fail_filters", "ccs_zmws_shortcut_filters",
                  "ccs_zmws_pass_filters_pct", "ccs_zmws_fail_filters_pct", "ccs_zmws_shortcut_filters_pct", "ref_map",
                  "gcs_input_dir", "subreads_bam", "subreads_pbi", "ccs_bam", "ccs_pbi", "hifi_bam", "hifi_pbi",
                  "input_bam", "input_pbi"]
    tbl_rows = []
    for e in tqdm(ts):
        r = load_ccs_report(project, e['Files']['ccs_reports.txt'], e)

        application = e['WellSample'][0]['Application'] if 'Application' in e['WellSample'][0] else "longReads"
        experiment_type = "CLR"

        if ('IsCCS' in e['WellSample'][0] and e['WellSample'][0]['IsCCS'] == 'true') or e['Files']['reads.bam'] != "" or application == 'hifiReads':
            experiment_type = "CCS"
        if ('IsoSeq' in e['WellSample'][0] and e['WellSample'][0]['IsoSeq'] == 'true') or (application == 'isoSeq'):
            experiment_type = "ISOSEQ"
        if ('BioSample' in e and 'MAS' in e['BioSample'][0]['Name']) or ('custom' in application):
            experiment_type = "MASSEQ"

        if e['Files']['hifi_reads.bam'] != "":
            input_bam = e['Files']['hifi_reads.bam']
            input_pbi = e['Files']['hifi_reads.bam.pbi']
        elif e['Files']['reads.bam'] != "":
            input_bam = e['Files']['reads.bam']
            input_pbi = e['Files']['reads.bam.pbi']
        elif e['Files']['subreads.bam'] != "":
            input_bam = e['Files']['subreads.bam']
            input_pbi = e['Files']['subreads.bam.pbi']

        #correctable = False if e['Files']['subreads.bam'] == "" else is_correctable(e['Files']['subreads.bam'])
        #if correctable:
        #    e['WellSample'][0]['IsCCS'] = 'true'

        tbl_rows.append([
            e['CollectionMetadata'][0]['UniqueId'] if 'Context' in e['CollectionMetadata'][0] else "",

            e['CellPac'][0]['Barcode'] if 'Barcode' in e['CellPac'][0] else "UnknownFlowcell",

            e['CollectionMetadata'][0]['InstrumentName'] if 'Context' in e['CollectionMetadata'][0] else "UnknownInstrument",
            e['CollectionMetadata'][0]['Context'] if 'Context' in e['CollectionMetadata'][0] else "UnknownMovie",
            e['WellSample'][0]['WellName'] if 'WellName' in e['WellSample'][0] else "Z00",
            e['WellSample'][0]['CreatedAt'] if 'CreatedAt' in e['WellSample'][0] else "0001-01-01T00:00:00",
            re.sub("[# ]", "", e['BioSample'][0]['Name']) if 'BioSample' in e else "UnknownBioSample",
            re.sub("[# ]", "", e['WellSample'][0]['Name']) if 'Name' in e['WellSample'][0] else "UnknownWellSample",
            e['WellSample'][0]['InsertSize'] if 'InsertSize' in e['WellSample'][0] else "0",
            e['WellSample'][0]['IsCCS'] if 'IsCCS' in e['WellSample'][0] else "unknown",
            e['WellSample'][0]['IsoSeq'] if 'IsoSeq' in e['WellSample'][0] else "unknown",

            "true" if 'ConsensusReadSet' in e else "false",
            e['WellSample'][0]['Description'] if 'Description' in e['WellSample'][0] else "unknown",
            e['WellSample'][0]['Application'] if 'Application' in e['WellSample'][0] else "unknown",
            experiment_type,

            e['DataSetMetadata'][0]['NumRecords'],
            e['DataSetMetadata'][0]['TotalLength'],

            e['Files']['ccs_reports.txt'],
            r['ZMWs input'],
            r['ZMWs pass filters'],
            r['ZMWs fail filters'],
            r['ZMWs shortcut filters'],
            "{:.2f}".format(100.0 * r['ZMWs pass filters'] / (
                    r['ZMWs pass filters'] + r['ZMWs fail filters'] + r['ZMWs shortcut filters'] + 1)),
            "{:.2f}".format(100.0 * r['ZMWs fail filters'] / (
                    r['ZMWs pass filters'] + r['ZMWs fail filters'] + r['ZMWs shortcut filters'] + 1)),
            "{:.2f}".format(100.0 * r['ZMWs shortcut filters'] / (
                    r['ZMWs pass filters'] + r['ZMWs fail filters'] + r['ZMWs shortcut filters'] + 1)),

            'gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/grch38_noalt.txt',

            e['Files']['input_dir'],
            e['Files']['subreads.bam'],
            e['Files']['subreads.bam.pbi'],
            e['Files']['reads.bam'],
            e['Files']['reads.bam.pbi'],
            e['Files']['hifi_reads.bam'],
            e['Files']['hifi_reads.bam.pbi'],

            input_bam,
            input_pbi
        ])

    tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)
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

    joined_tbl['description'] = joined_tbl['description'].str.replace(r'\s+', ' ', regex=True).astype('str')
    joined_tbl['bio_sample'] = joined_tbl['bio_sample'].str.replace(r'\s+', ' ', regex=True).astype('str')
    joined_tbl['well_sample'] = joined_tbl['well_sample'].str.replace(r'\s+', ' ', regex=True).astype('str')

    return joined_tbl


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


def load_xmls(gcs_buckets, project):
    storage_client = storage.Client(project=project)

    input_dirs = set()
    for gcs_bucket in gcs_buckets:
        bucket_name = re.sub("^gs://", "", gcs_bucket)
        bucket = storage_client.bucket(bucket_name)

        blobs = storage_client.list_blobs(bucket)
        for blob in blobs:
            if '.subreadset.xml' in blob.name or \
               '.consensusreadset.xml' in blob.name or \
               '.subreads.bam' in blob.name or \
               '.reads.bam' in blob.name or \
               '.hifi_reads.bam' in blob.name:

                input_dirs.add((bucket_name, os.path.dirname(blob.name)))

    ts = []
    for bucket_name, input_dir in input_dirs:
        xml_file = None

        blobs = storage_client.list_blobs(bucket_name, prefix=input_dir)
        for blob in blobs:
            if '.consensusreadset.xml' in blob.name or '.subreadset.xml' in blob.name:
                xml_file = blob.name
                break

        if xml_file is not None:
            xml_blob = storage.Blob(bucket=bucket, name=xml_file)
            xml = xml_blob.download_as_string()
            doc = xmltodict.parse(xml)

            t = combine(traverse_xml('root', doc))
            t['Files'] = {
                'input_dir': input_dir,

                'subreadset.xml': "",
                'consensusreadset.xml': "",

                'subreads.bam': "",
                'subreads.bam.pbi': "",
                'reads.bam': "",
                'reads.bam.pbi': "",
                'hifi_reads.bam': "",
                'hifi_reads.bam.pbi': "",

                'ccs_reports.txt': ""
            }

            blobs = storage_client.list_blobs(bucket_name, prefix=input_dir)
            for blob in blobs:
                if '.subreadset.xml' in blob.name:
                    t['Files']['subreadset.xml'] = f'gs://{bucket_name}/{blob.name}'

                elif '.consensusreadset.xml' in blob.name:
                    t['Files']['consensusreadset.xml'] = f'gs://{bucket_name}/{blob.name}'

                elif '.subreads.bam.pbi' in blob.name:
                    t['Files']['subreads.bam.pbi'] = f'gs://{bucket_name}/{blob.name}'

                elif '.subreads.bam' in blob.name and '.pbi' not in blob.name:
                    t['Files']['subreads.bam'] = f'gs://{bucket_name}/{blob.name}'

                elif '.reads.bam.pbi' in blob.name:
                    t['Files']['reads.bam.pbi'] = f'gs://{bucket_name}/{blob.name}'

                elif '.reads.bam' in blob.name and '.pbi' not in blob.name:
                    t['Files']['reads.bam'] = f'gs://{bucket_name}/{blob.name}'

                elif '.hifi_reads.bam.pbi' in blob.name:
                    t['Files']['hifi_reads.bam.pbi'] = f'gs://{bucket_name}/{blob.name}'

                elif '.hifi_reads.bam' in blob.name and '.pbi' not in blob.name:
                    t['Files']['hifi_reads.bam'] = f'gs://{bucket_name}/{blob.name}'

                elif '.ccs_reports.txt' in blob.name:
                    t['Files']['ccs_reports.txt'] = f'gs://{bucket_name}/{blob.name}'

            ts.append(t)

    return ts


def load_ccs_report(project, ccs_report_path, e):
    d = {
        'ZMWs input': 0,
        'ZMWs pass filters': 0,
        'ZMWs fail filters': 0,
        'ZMWs shortcut filters': 0,
        'ZMWs with tandem repeats': 0,
        'Below SNR threshold': 0,
        'Median length filter': 0,
        'Lacking full passes': 0,
        'Heteroduplex insertions': 0,
        'Coverage drops': 0,
        'Insufficient draft cov': 0,
        'Draft too different': 0,
        'Draft generation error': 0,
        'Draft above --max-length': 0,
        'Draft below --min-length': 0,
        'Reads failed polishing': 0,
        'Empty coverage windows': 0,
        'CCS did not converge': 0,
        'CCS below minimum RQ': 0,
        'Unknown error': 0
    }

    if ccs_report_path != "":
        storage_client = storage.Client(project=project)

        ccs_report = re.sub("^gs://", "", e['Files']['ccs_reports.txt']).split("/")
        blobs = storage_client.list_blobs(ccs_report[0], prefix="/".join(ccs_report[1:]))

        for blob in blobs:
            blob.download_to_filename("ccs_report.txt")

            f = open("ccs_report.txt", "r")

            d = {}
            for line in f:
                if len(line) > 1 and 'Exclusive' not in line and 'Additional' not in line:
                    a = line.rstrip().split(":")

                    try:
                        k = a[0].rstrip()
                        v = float(re.sub(" ", "", re.sub(" \(.*$", "", a[1])))

                        if k not in d:
                            d[k] = 0.0;

                        d[k] = d[k] + v
                    except:
                        pass

            break

    return d


def is_correctable(subreads_bam):
    cmd = 'gsutil cat ' + subreads_bam + ' | samtools view | head -n 1000 | cut -d"/" -f2 | uniq -c'
    lines = subprocess.getoutput(cmd).split("\n")

    num_ccs_reads = 0
    num_reads = 1 # add a pseudocount to avoid div by 0
    for l in lines:
        l = l.strip()

        if l != "":
            a = l.strip().split(" ")
            try:
                if int(a[0]) >= 3:
                    num_ccs_reads += 1
                num_reads += 1
            except:
                pass

    return num_ccs_reads / num_reads >= 0.3


def update_sample_table(namespace, workspace, buckets, project, tbl_old):
    tbl_new = load_new_sample_table(buckets, project)
    joined_tbl = merge_tables(tbl_old, tbl_new)
    joined_tbl = joined_tbl.replace('^nan$', '', regex=True)

    return joined_tbl


def update_sample_set_table(namespace, workspace, ss_old, membership, joined_tbl):
    # create old membership set
    oms = pd \
        .DataFrame({'entity:sample_set_id': list(ss_old['entity:sample_set_id']), 'sample': membership}) \
        .explode('sample', ignore_index=True)
    oms.columns = ['membership:sample_set_id', 'sample']

    # create sample set
    ss = joined_tbl.filter(['bio_sample'], axis=1).drop_duplicates()
    ss.columns = [f'entity:sample_set_id']

    if ss_old is not None:
        ss = pd.merge(ss_old, ss, how='outer', sort=True)
    ss = ss.replace('^nan$', '', regex=True)

    # create new membership set
    ms = joined_tbl.filter(['bio_sample', 'entity:sample_id'], axis=1).drop_duplicates()
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
    parser = argparse.ArgumentParser(description='Update Terra workspace sample table', prog='update_pacbio_tables')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-r', '--run', action='store_true', help="Turn off the default dry-run mode")
    parser.add_argument('buckets', metavar='B', type=str, nargs='+', help='GCS buckets to scan')
    args = parser.parse_args()

    tbl_old, _ = load_table(args.namespace, args.workspace, 'sample')
    ss_old, membership = load_table(args.namespace, args.workspace, 'sample_set', store_membership=True)

    s = update_sample_table(args.namespace, args.workspace, args.buckets, args.project, tbl_old)
    ss, nms = update_sample_set_table(args.namespace, args.workspace, ss_old, membership, s)

    if args.run:
        fapi.__SESSION = None
        upload_tables(args.namespace, args.workspace, s, ss, nms)


if __name__ == "__main__":
    main()
