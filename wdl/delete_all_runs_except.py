#!/usr/bin/env python3

import os

os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
from subprocess import check_output


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workflow_name")
    p.add_argument("workflow_id", nargs="+")
    args = p.parse_args()

    workflows_path = f"gs://broad-dsde-methods/cromwell-execution-34/{args.workflow_name}"
    directories = check_output(f"gsutil ls {workflows_path}", shell=True)
    directories = directories.decode('UTF-8').strip()

    if not directories:
        p.error(f"No subdirs found in {workflows_path}")

    directories = [d.strip('/') for d in directories.split("\n")]
    print(f"==> Found {len(directories)} subdirs in {workflows_path}")
    existing_workflow_ids = {os.path.basename(d) for d in directories}
    workflows_ids_to_keep = set(args.workflow_id)

    unknown_ids = workflows_ids_to_keep - existing_workflow_ids
    if len(unknown_ids) > 0:
        print("\n".join(f"{workflows_path}/{wid}" for wid in unknown_ids))
        p.error(f"Directories don't exist for these ids: {', '.join(unknown_ids)}")

    workflows_ids_to_delete = existing_workflow_ids - workflows_ids_to_keep
    if not workflows_ids_to_delete:
        p.exit("No additional directories to delete. Exiting..\n")

    print(f"Will keep these {len(workflows_ids_to_keep)} directories:")
    print("\n".join(f"{workflows_path}/{wid}" for wid in workflows_ids_to_keep))

    print(f"Will delete these {len(workflows_ids_to_delete)} directories:")
    print("\n".join(f"{workflows_path}/{wid}" for wid in workflows_ids_to_delete))
    if input("Continue? [Y/n] ").upper() != "Y":
        p.exit(message="Exiting..\n")

    os.system("gsutil -m rm -rf " + " ".join(f"{workflows_path}/{wid}" for wid in workflows_ids_to_delete))
    print("Done")
