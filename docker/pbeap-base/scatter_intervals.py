#!/usr/bin/env python3

import argparse
import os
import re

p = argparse.ArgumentParser(description="Takes a .bed file path and splits it into the specified number of subsets")

g = p.add_mutually_exclusive_group(required=True)
g.add_argument("--num-records-per-shard", type=int, help="Write at most this many records to each output .bed file")
g.add_argument("--num-shards", help="Create this many output .bed files", type=int)

p.add_argument("--interleaved", action="store_true", help="Split input records into shards in an interleaved pattern "
    "(123123123123 instead of 111122223333 for 3 shards).")
p.add_argument("--output-dir-suffix", default="_sharded", help="Output .bed files will be written to an output "
    " directory that's named like the input .bed file name (without '.bed') + this suffix.")
p.add_argument("bed_file_path", nargs="+", help="One or more input .bed files to split. Each will be split separately.")

args = p.parse_args()

for bed_file_path in args.bed_file_path:
    # compute num_records in bed_file_path
    with open(bed_file_path) as f:
        num_records = sum(1 for line in f)

    # compute args.num_shards or args.num_records_per_shard when the other is specified
    if args.num_records_per_shard:
        args.num_shards = (num_records // args.num_records_per_shard) + 1
    elif args.num_shards:
        args.num_records_per_shard = (num_records // args.num_shards) + 1

    # create output dir and delete previous files
    basename = os.path.basename(bed_file_path)
    output_filename_prefix = re.sub(".bed$", "", basename)
    output_dir = f"{output_filename_prefix}{args.output_dir_suffix}"
    os.system(f"mkdir -p {output_dir}")
    os.system(f"rm {output_dir}/{output_filename_prefix}.*.bed")

    # split input bed_file into shards
    output_files = {}
    print(
        f"==> {bed_file_path} has {num_records} records. Generating {args.num_shards} files in {output_dir}, "
        f"with {args.num_records_per_shard} records in each.")

    with open(bed_file_path) as f:
        for i, line in enumerate(f):
            if args.interleaved:
                shard_i = i % args.num_shards
            else:
                shard_i = i // args.num_records_per_shard

            output_filename = f"{output_filename_prefix}.{shard_i:02d}.bed"

            if output_filename not in output_files:
                print(f"Creating {output_dir}/{output_filename}")
                output_files[output_filename] = open(f"{output_dir}/{output_filename}", "w")

            output_files[output_filename].write(line)

    # close shard files
    for output_file in output_files.values():
        output_file.close()

    print("Done")
