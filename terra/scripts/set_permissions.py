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


def main():
    parser = argparse.ArgumentParser(description='Set permissions for a user to a workspace', prog='set_permissions')
    parser.add_argument('-n', '--namespace', type=str, help="Terra namespace")
    parser.add_argument('-w', '--workspace', type=str, help="Terra workspace")
    parser.add_argument('-u', '--user', type=str, action="append", help="Usernames to add as owners to the workspace")
    parser.add_argument('-l', '--level', type=str, help="Permission level to set (e.g. READER, WRITER, OWNER)")
    args = parser.parse_args()

    print(f"Accessing Terra as '{fapi.whoami()}'.")

    owners = [{"email": u, "accessLevel": args.level} for u in args.user]

    b = fapi.update_workspace_acl(args.namespace, args.workspace, owners)


if __name__ == "__main__":
    main()
