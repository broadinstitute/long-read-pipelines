#!/usr/bin/env python3

import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

import argparse
import logging
import sys

from run_cromwell_job_utils import execute

logging.basicConfig(level=logging.INFO)

FILE_PATH_PREFIX = {
    "local": "/Users/kiran/repositories/PBEAP/",
    "remote": "gs://broad-dsde-methods-kiran/",
}

WORKFLOWS = {
    "basic_stats": { "name": "basic_stats", "wdl": "/Users/kiran/repositories/PBEAP/wdl/basic_stats/basic_stats.wdl" },
}


REFERENCE_GENOMES = {
    "fasta": {
        "37": "resources/references/grch37/hs37d5.fa",
    },
}

PBEAP_SUBREAD_BAMS = {
    "37": {
        "Ecoli_Training_Run": { "path": "pb_eap/Ecoli_Training_Run/m64020_190123_195048.subreads.bam" }
    },
}

PBEAP_DATA = {
    bam_label: {
        "37": {
            "input_bam": bam_info["path"],
        }
    } for bam_label, bam_info in PBEAP_SUBREAD_BAMS["37"].items()
}

DATA = {}
DATA.update(PBEAP_DATA)

#LOCAL_WORKFLOW_OPTIONS = {}
GCS_WORKFLOW_OPTIONS = {
    "write_to_cache": True,
    "read_from_cache": True
}

MAX_SHARDS_LIMIT = {
    "local": 1,
    "remote": 1,

    "basic_stats": 1,
}

def parse_args():
    p = argparse.ArgumentParser()

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--local", action="store_true")
    group.add_argument("--remote", action="store_true")

    p.add_argument("-g", "--genome-version", choices=REFERENCE_GENOMES["fasta"].keys(), required=True)
    p.add_argument("--skip-input-validation", action="store_true", help="Skip input validation.")
    p.add_argument("--test", action="store_true", help="Don't actually submit the job.")

    p.add_argument("tool", choices=WORKFLOWS.keys())
    p.add_argument("data", nargs="+", choices=DATA.keys())
    #p.add_argument("label", nargs="?", help="optional label describing this run")

    args = p.parse_args()

    return args


def main():
    args = parse_args()

    tool_name = WORKFLOWS[args.tool]['name']
    wdl_path = WORKFLOWS[args.tool]['wdl']
    workflow_options = GCS_WORKFLOW_OPTIONS
    data_root_path = FILE_PATH_PREFIX["local"] if args.local else FILE_PATH_PREFIX["remote"]
    max_shards = min(
        MAX_SHARDS_LIMIT[args.tool],
        MAX_SHARDS_LIMIT["local" if args.local else "remote"],
    )


    # generate input json
    for data in args.data:
        if args.genome_version not in DATA[data]:
            sys.exit(f"{data} data not found for GRCh{args.genome_version}.")

        input_json = DATA[data][args.genome_version]
        input_json.update({
            "output_prefix": f"results/{args.tool}/grch{args.genome_version}/{data}",
            "ref_fasta": REFERENCE_GENOMES["fasta"][args.genome_version],
            "ref_fasta_fai": REFERENCE_GENOMES["fasta"][args.genome_version] + ".fai",
            "ref_fasta_dict": REFERENCE_GENOMES["fasta"][args.genome_version].replace(".fa", ".dict"),
        })

        if args.tool == "stretch":
            input_json["reference_files_dir"] = REFERENCE_GENOMES["stretch_ref_dir"][args.genome_version]

        elif args.tool in ("gatk", "gatk-str"):
            input_json["use_str_model"] = "true" if args.tool == "gatk-str" else "false"

        for key, value in input_json.items():
            if args.local and key == "num_shards":
                # limit num_shards to 3 for local
                input_json[key] = min(value, max_shards)


        execute(args, tool_name, f"hg{args.genome_version}_{data}", wdl_path, workflow_options, input_json, data_root_path)


if __name__ == "__main__":
    main()
