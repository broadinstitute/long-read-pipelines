#!/usr/bin/env python

import re
import sys
import os
import pysam

from tqdm import tqdm

if len(sys.argv) != 3:
    print(f"USAGE: {os.path.basename(__file__)} INPUT_BAM OUTPUT_BAM", file=sys.stderr)
    sys.exit(1)

bam_file_name = sys.argv[1]
out_bam_file_name = sys.argv[2]

if not os.path.exists(bam_file_name):
    print(f"ERROR: Given input file does not exist: {bam_file_name}", file=sys.stderr)
    sys.exit(1)

if os.path.exists(out_bam_file_name):
    print(f"ERROR: Given output file already exists: {out_bam_file_name}", file=sys.stderr)
    sys.exit(1)

print(f"Replacing tags in {bam_file_name}.")
print(f"    output file: {out_bam_file_name}.")

tag_replacement_map = {
    "Jp": "CB",
    "Jq": "ZU",
    "Jr": "XM",
    "Js": "XC",
    "Jt": "it",
}

print()
print("Tag replacement map:")
for old, new in tag_replacement_map.items():
    print(f"    {old} -> {new}")
print()

with pysam.AlignmentFile(bam_file_name, "rb", check_sq=False, require_index=False) as bam_file:
    with pysam.AlignmentFile(out_bam_file_name, "wb", header=bam_file.header) as out_bam_file:
        for i, read in enumerate(tqdm(bam_file, desc="fixing tags", unit=" read")):
            for old, new in tag_replacement_map.items():
                if read.has_tag(old):
                    read.set_tag(new, read.get_tag(old))
                    read.set_tag(old, None)

                out_bam_file.write(read)
