#!/usr/bin/env python

import sys
import time
import re
import pathlib

import pysam

from tqdm import tqdm

ZMW_TAG = 'zm'
zmw_regex = re.compile(r'.*?/(\d+)/.*')

if __name__ == '__main__':

    # Right now we accept only 1 argument - the directory in which the quant files live:
    if len(sys.argv) != 2:
        print(f"{sys.argv[0]} MAS_SEQ_ARRAY_ELEMENT_PB_BAM", file=sys.stderr)
        print(f"Extracts one read from each ZMW from the given MAS-seq array element bam file.", file=sys.stderr)
        print(f"Assumes reads originated from a PacBio instrument and a MAS-seq library, and that they have already been split into array elements.", file=sys.stderr)
        print(f"Also Assumes that the reads are sorted by read name.", file=sys.stderr)
        print(f"", file=sys.stderr)
        print(f"Error: you must specify a pac bio bam file to downsample", file=sys.stderr)
        sys.exit(1)

    in_file_path = pathlib.Path(sys.argv[1])

    if not in_file_path.exists():
        print(f"ERROR: File does not exist: {in_file_path}", file=sys.stderr)
        sys.exit(1)
    if not in_file_path.is_file():
        print(f"ERROR: Given input file is not a file.  Is it a directory?", file=sys.stderr)
        sys.exit(1)
    if not (in_file_path.suffix == ".sam" or in_file_path.suffix == ".bam"):
        print(f"ERROR: Only bam/sam files are supported.", file=sys.stderr)
        sys.exit(1)

    is_bam = in_file_path.suffix == '.bam'

    out_file_name = in_file_path.stem + ".ZMW_downsampled.bam"

    last_timing = time.time()
    read_count = 0

    zmw_set = set()

    has_XN_tag = False
    is_first = True

    num_written = 0
    num_skipped = 0

    with pysam.AlignmentFile(in_file_path, 'rb' if is_bam else 'r', check_sq=False) as bam_file:
        with pysam.AlignmentFile(out_file_name, 'wb', check_sq=False, header=bam_file.header) as output_file:
            with tqdm(desc=f"Downsampling bam file", unit=" read") as pbar:
                for read in bam_file.fetch(until_eof=True):

                    if is_first and read.has_tag("XN"):
                        print("Read has XN Tag.  Using XN for names!")
                        has_XN_tag = True
                    is_first = False

                    name = read.get_tag("XN") if has_XN_tag else read.query_name

                    i1 = name.find("/")
                    i2 = name.find("/", i1+1)
                    zmw = int(name[i1+1:i2])

                    if zmw not in zmw_set:
                        output_file.write(read)
                        zmw_set.add(zmw)
                        num_written += 1
                    else:
                        num_skipped += 1

                    read_count += 1
                    pbar.update(1)

print(f"Total reads: {read_count}")
print(f"Num reads output: {num_written}")
print(f"Num reads skipped: {num_skipped}")


