#!/usr/bin/env python

import os
import re
import math
import hashlib
import argparse
import copy

import numpy as np
import pandas as pd
import firecloud.api as fapi

from google.cloud import storage

from collections import OrderedDict

from tqdm import tqdm
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


def load_summaries(gcs_buckets):
    storage_client = storage.Client()
    schemas = OrderedDict()

    ts = {}
    for gcs_bucket in gcs_buckets:
        blobs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket), prefix="inputs")

        for blob in blobs:
            if blob.name.endswith(".fast5") and "fail" not in blob.name:
                gcs_path = os.path.dirname(os.path.dirname(blob.name))
                
                ts[gcs_path] = blob.time_created

    return ts


def load_summaries_old(gcs_buckets, project):
    storage_client = storage.Client(project=project)

    ts = []
    for gcs_bucket in gcs_buckets:
        blobs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket), timeout=120)

        for blob in blobs:
            if 'final_summary' in blob.name:
                doc = blob.download_as_string(timeout=120)
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
                                               prefix=os.path.dirname(blob.name) + "/" + t['sequencing_summary_file'],
                                               timeout=120)
                for b in bs:
                    t['Files']['sequencing_summary.txt'] = gcs_bucket + "/" + b.name

                if 'sequencing_summary.txt' not in t['Files']:
                    pp = pprint.PrettyPrinter(indent=4)
                    pp.pprint(t)

                # Handle barcoded datasets
                cs = storage_client.list_blobs(re.sub("^gs://", "", gcs_bucket),
                                               prefix=re.sub("^/", "", re.sub(gcs_bucket, "", t['Files']['fastq_pass_dir'])),
                                               timeout=120)
                dirs = set()
                for c in cs:
                    # Barcode identifier
                    d = os.path.basename(os.path.dirname(c.name))

                    if d != "fastq_pass":
                        dirs.add(d)

                for d in dirs:
                    tb = copy.deepcopy(t)
                    tb['Files']['fastq_pass_dir'] = t['Files']['fastq_pass_dir'] + "/" + d
                    tb['sample_id'] = t['sample_id'] + "." + d
                    tb["entity:sample_id"] = hashlib.md5(tb["Files"]["fastq_pass_dir"].encode("utf-8")).hexdigest()

                    ts.append(tb)

                t["entity:sample_id"] = hashlib.md5(t["Files"]["final_summary.txt"].encode("utf-8")).hexdigest()
                ts.append(t)

    return ts


def load_new_sample_table(default_bucket, project):
    configs = {
        'FLO-MIN106': {'SQK-LSK109': 'dna_r9.4.1_450bps_sup.cfg'},
        'FLO-MIN111': {'SQK-LSK109': 'dna_r9.4.1_450bps_sup.cfg',
                       'SQK-LSK110': 'dna_r10.3_450bps_sup.cfg'},
        'FLO-MIN112': {'SQK-LSK112': 'dna_r10.4_e8.1_sup.cfg'},
        'FLO-PRO002': {'SQK-LSK109': 'dna_r9.4.1_450bps_sup_prom.cfg'}
    }

    ts = load_summaries([default_bucket])

    storage_client = storage.Client()
    
    tbl_rows = []
    for k,v in ts.items():
        if k != "inputs":
            f = default_bucket + "/" + k
    
            fs, ss = "unknown", "unknown"
            d = {
                'instrument': 'unknown', #instrument=MC-110675
                'position': 'unknown', #position=MC-110675_0
                'protocol_group_id': 'unknown', #protocol_group_id=coi2_17may2021
                'flow_cell_id': 'unknown', #flow_cell_id=FAO99587
                'sample_id': 'unknown', #=no_sample
                'protocol': 'unknown', #=sequencing/sequencing_MIN106_DNA:FLO-MIN106:SQK-LSK109
                'basecalling_enabled': 'unknown', #==1
            }
    
            blobs = storage_client.list_blobs(re.sub("^gs://", "", default_bucket), prefix=k)
            for b in blobs:
                if re.search("\/final_summary.*.txt", b.name):
                    fs = default_bucket + "/" + b.name
    
                    q = b.download_as_text().splitlines()
    
                    for a in q:
                        if '=' in a:
                            k2,v2 = a.split("=")
    
                            if v2 != "":
                                d[k2] = v2
    
                if re.search("\/sequencing_summary.*.txt", b.name):
                    ss = default_bucket + "/" + b.name
    
            basecalling_model = ''
            a = re.sub('sequencing/sequencing_', '', d['protocol']).split(":")
            if len(a) == 3:
                inst_type, flowcell_type, kit_type = a
    
                basecalling_model = configs[flowcell_type][kit_type]
                    
            tbl_rows.append([
                hashlib.md5(f.encode("utf-8")).hexdigest(),
                f,
                fs,
                ss,
                d['instrument'],
                d['position'],
                d['protocol_group_id'],
                d['flow_cell_id'],
                d['protocol'],
                str(v),
                d['sample_id'],
                basecalling_model,
                "",
                "none",
                os.path.basename(f)
            ])
        
    tbl_header = ["entity:sample_id", "fast5_dir", "final_summary", "sequencing_summary", "instrument", "position",
                  "protocol_group_id", "flow_cell_id", "protocol", "upload_date", "sample_name", "basecalling_model",
                  "barcode_kit", "notes", "fast5_dir_basename"]
        
    tbl_new = pd.DataFrame(tbl_rows, columns=tbl_header)
    tbl_new = tbl_new.sort_values(['fast5_dir_basename', 'upload_date'])
    tbl_new = tbl_new.groupby('fast5_dir_basename').first().reset_index().drop(columns='fast5_dir_basename')
    
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

    joined_tbl = merge_tables(tbl_old, tbl_new) if tbl_old is not None else tbl_new
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
    ss = ss.replace('^nan$', '', regex=True)

    # create new membership set
    ms = joined_tbl.filter(['sample_name', 'entity:sample_id'], axis=1).drop_duplicates()
    ms.columns = [f'membership:sample_set_id', f'sample']

    # create full membership set
    fms = pd.merge(ms, oms, how='outer', indicator=True)

    # create new/modified membership set
    nms = fms[fms['_merge'] != 'both']

    return ss, nms


def upload_table(namespace, workspace, sample_tbl):
    columns_reordered = ['entity:sample_id'] + list(filter(lambda x: x != 'entity:sample_id', list(sample_tbl.columns)))
    sample_tbl = sample_tbl[columns_reordered]
    
    a = fapi.upload_entities(namespace, workspace, entity_data=sample_tbl.to_csv(index=False, sep="\t"), model='flexible')
    
    if a.status_code == 200:
        print(f'Uploaded {len(sample_tbl)} rows successfully.')
    else:
        print(a.json())


def main():
    parser = argparse.ArgumentParser(description='Update Terra cigass-oxfordnano sample table', prog='update_cigass_oxfordnano_tables')
    parser.add_argument('-p', '--project', type=str, help="GCP project")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-r', '--run', action='store_true', help="Turn off the default dry-run mode")
    parser.add_argument('bucket', metavar='B', type=str, help='GCS bucket to scan')
    args = parser.parse_args()

    s = update_sample_table(args.namespace, args.workspace, args.bucket, args.project)

    if args.run:
        upload_table(args.namespace, args.workspace, s)


if __name__ == "__main__":
    main()
