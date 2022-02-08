#!/usr/bin/env python

import time
import sys
import pysam

# Right now we require 2 args:
if len(sys.argv) != 3:
    print(f"{sys.argv[0]} BAM_FILE OUT_BASE_NAME", file=sys.stderr)
    print(f"Error: you must specify an input BAM_FILE and an OUT_BASE_NAME.", file=sys.stderr)
    sys.exit(1)

pysam.set_verbosity(0)

in_bam_file = sys.argv[1]
out_bam_base_name = sys.argv[2]

print(f"Processing file: {in_bam_file}", file=sys.stderr)

start_t = time.time()

with pysam.AlignmentFile(in_bam_file, "rb", check_sq=False, require_index=False) as bam_file:
    last_contig_name = None
    out_bam = None
    for read in bam_file:
        if read.reference_name != last_contig_name:
            print(f"New Contig: {read.reference_name} ({out_bam_base_name}.{read.reference_name}.bam)", file=sys.stderr)
            if out_bam:
                out_bam.close()
            out_bam = pysam.AlignmentFile(f"{out_bam_base_name}.{read.reference_name}.bam", "wb", header=bam_file.header)
            last_contig_name = read.reference_name
        out_bam.write(read)
    if out_bam:
        out_bam.close()

print('Done!', file=sys.stderr)
end_t = time.time()
print(f"Elapsed time: {end_t - start_t}s", file=sys.stderr)
