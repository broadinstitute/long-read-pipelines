#!/usr/bin/env python

import sys
import time
import re
import pathlib

import pysam

ZMW_TAG = 'zm'
zmw_regex = re.compile(r'.*?/(\d+)/.*')

if __name__ == '__main__':

    # Right now we accept only 1 argument - the directory in which the quant files live:
    if len(sys.argv) != 2:
        print(f"{sys.argv[0]} PAC_BIO_BAM", file=sys.stderr)
        print(f"Adds the ZMW tag to each read in the given bam file if it does not already exist.", file=sys.stderr)
        print(f"Assumes read names conform to the pacbio standard readname (MOVIE/ZMW/[ccs|molecular coordinates])", file=sys.stderr)
        print(f"", file=sys.stderr)
        print(f"Error: you must specify a pac bio bam file to update.", file=sys.stderr)
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

    out_file_name = in_file_path.stem + ".with_zmws.bam"

    last_timing = time.time()
    read_count = 0
    with pysam.AlignmentFile(in_file_path, 'rb' if is_bam else 'r', check_sq=False) as bam_file:
        with pysam.AlignmentFile(out_file_name, 'wb', check_sq=False, header=bam_file.header) as output_file:
            for read in bam_file.fetch(until_eof=True):
                if read_count % 1000 == 0:
                    current_time = time.time()
                    print('Processed reads: {:,}. Time elapsed: {:.2f}s'.format(read_count, current_time - last_timing), file=sys.stderr)
                    last_timing = current_time

                try:
                    zmw = int(read.get_tag(ZMW_TAG))
                except KeyError:
                    try:
                        zmw = int(zmw_regex.match(read.query_name).group(1))
                    except AttributeError:
                        print(f"FATAL ERROR: Could not find ZMW in read {read_count}: {read.query_name}",
                              file=sys.stderr)
                        sys.exit(2)

                read.set_tag(ZMW_TAG, zmw, value_type='i')
                output_file.write(read)

                read_count += 1
