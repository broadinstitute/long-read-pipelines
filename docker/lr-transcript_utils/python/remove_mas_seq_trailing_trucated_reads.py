#!/usr/bin/env python

import sys
import time
import re
import pathlib

import pysam

READ_ALTERED_NAME_TAG = "XN"
zmw_regex = re.compile(r'.*?/(\d+)/.*')

if __name__ == '__main__':

    # Right now we accept exactly 2 arguments:
    if len(sys.argv) != 3:
        print(f"{sys.argv[0]} LONGBOW_SEGMENTED_BAM PREFIX", file=sys.stderr)
        print(f"Removes any array element with the -END suffix, indicating that it is a trailing truncated read.", file=sys.stderr)
        print(f"", file=sys.stderr)
        print(f"Error: you must specify a LONGBOW_SEGMENTED_BAM file to update and a PREFIX for the output file.", file=sys.stderr)
        sys.exit(1)

    in_file_path = pathlib.Path(sys.argv[1])
    prefix = sys.argv[2]

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

    out_file_name = prefix + ".bam"

    last_timing = time.time()
    read_count = 0
    num_end_elements_removed = 0
    with pysam.AlignmentFile(in_file_path, 'rb' if is_bam else 'r', check_sq=False) as bam_file:
        with pysam.AlignmentFile(out_file_name, 'wb', check_sq=False, header=bam_file.header) as output_file:
            for read in bam_file.fetch(until_eof=True):
                if read_count % 1000 == 0:
                    current_time = time.time()
                    print('Processed reads: {:,}. Time elapsed: {:.2f}s'.format(read_count, current_time - last_timing), file=sys.stderr)
                    last_timing = current_time

                # Make sure we're looking at the correct name.
                # After reads are segmented, the longbow-specific name gets put into the READ_ALTERED_NAME_TAG tag.
                # As a fall-back, we can check the read name itself as well:
                name = read.get_tag(READ_ALTERED_NAME_TAG) if read.has_tag(READ_ALTERED_NAME_TAG) else read.query_name

                if name.endswith("-END"):
                    num_end_elements_removed += 1
                    read_count += 1
                    continue
                else:
                    output_file.write(read)

                read_count += 1

    print(f"Num reads: {read_count}")
    print(f"Num end element reads removed: {num_end_elements_removed} ({100*num_end_elements_removed/read_count:2.4f}%)")
