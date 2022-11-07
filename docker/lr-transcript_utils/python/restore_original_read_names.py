#!/usr/bin/env scripts

# This is a quick script to swap the read name with the value of the XN tag.
#
# Author: Jonn Smith
# Date:  2022 02 19

import sys
import time
import pysam

NAME_TAG = "XN"

if len(sys.argv) != 2:
    print("ERROR: Must provide only 1 argument - an input bam file.", file=sys.stderr)
    sys.exit(1)

# OK - read in the file
t_start = time.time()
input_bam = sys.argv[1]
print("Converting read names back to original names...")
with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:
    with pysam.AlignmentFile(f"out.original_read_names_restored.bam", "wb", header=bam_file.header) as out_bam:
        for i, read in enumerate(bam_file):

            if i % 10000 == 0:
                print(f"Processed {i} reads...")

            original_name = read.get_tag(NAME_TAG)
            hashed_name = read.query_name

            read.query_name = original_name
            read.set_tag(NAME_TAG, hashed_name)

            out_bam.write(read)

t_end = time.time()
print(f"Done!  Elapsed time: {t_end - t_start}s")
