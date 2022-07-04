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


def main():
    parser = argparse.ArgumentParser(description='Backup Terra workspace table(s)', prog='backup_workspace_tables')
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-o', '--outdir', type=str, help="Backup directory")
    args = parser.parse_args()

    tbl_old, _ = load_table(args.namespace, args.workspace, 'sample')

    path = f"{args.outdir}/{args.namespace}/{args.workspace}/"
    os.makedirs(path, mode = 0o777, exist_ok = True)

    tbl_old.to_csv(f"{path}/sample.csv", index=False)

    os.system(f"cd {path}; git add sample.csv; git commit -m 'Backup.'")


if __name__ == "__main__":
    main()
