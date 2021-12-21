#!/usr/bin/env python

import argparse
import re
import json
import requests

from os.path import exists

import pandas as pd
import firecloud.api as fapi


def main():
    parser = argparse.ArgumentParser(description='Send Slack notification on new data', prog='notify_on_new_data')
    parser.add_argument('-u', '--webhook-url', type=str, help="Slack webhook")
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-d', '--table-dir', type=str, help="Directory containing old tables")
    args = parser.parse_args()

    tbl_old, _ = load_table(args.namespace, args.workspace, 'sample')
    tbl_filtered = tbl_old[~tbl_old.workspace.isin(['nan', ''])]

    for rw in tbl_filtered['workspace'].unique():
        qa = fapi.get_workspace(args.namespace, rw)
        if qa.status_code == 200:
            tb, _ = load_table(args.namespace, rw, 'sample')
            ss, _ = load_table(args.namespace, rw, 'sample_set', store_membership=False)

            new_fcs = diff(tb, rw, 'sample', args.table_dir)
            new_sms = diff(ss, rw, 'sample_set', args.table_dir)

            if len(new_fcs) > 0 or len(new_sms) > 0:
                slack_data = {'text': f"New data delivered to {rw}: {len(new_fcs)} flowcells, {len(new_sms)} samples ({', '.join(new_sms)})."}
                print(slack_data)

                response = requests.post(
                    url = args.webhook_url,
                    data = json.dumps(slack_data),
                    headers = {'Content-Type': 'application/json'}
                )

                if response.status_code != 200:
                    raise ValueError(
                        'Request to slack returned an error %s, the response is:\n%s'
                        % (response.status_code, response.text)
                )


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


def diff(tbl, workspace, name, dir):
    cols = [v for v in tbl.columns if 'entity:' in v]

    old_tbl_path = re.sub("\s+", "_", f'{dir}/{workspace}.{name}.csv')

    old_tbl = []
    if exists(old_tbl_path):
        old_tbl = pd.read_csv(old_tbl_path)
        old_tbl = list(old_tbl[cols[0]])

    new_tbl = list(tbl[cols[0]])
    de = list(set(new_tbl) - set(old_tbl))

    tbl.to_csv(old_tbl_path, index=False)

    return de


if __name__ == "__main__":
    main()