#!/usr/bin/env python3

import argparse
import collections
import os
import sys

p = argparse.ArgumentParser(description="Generates a .bed file that can be passed to PrintReads.")
p.add_argument("--combine-intervals", action="store_true", help="Combine all intervals in each chromosome into 1 "
    "interval, which starts at the beginning of the 1st interval in that chromosome, and ends at the end of the last interval.")
p.add_argument("bed_file_path", help="Input .bed file path")

args = p.parse_args()

if not os.path.isfile(args.bed_file_path):
    p.exit(f"File not found: {args.bed_file_path}")

records = []
with open(args.bed_file_path) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.strip().split()
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        records.append((chrom, start, end))

if args.combine_intervals:
    min_coord = collections.defaultdict(int)
    max_coord = collections.defaultdict(int)
    for chrom, start, end in records:
        if not min_coord[chrom] or start < min_coord[chrom]:
            min_coord[chrom] = start
        if not max_coord[chrom] or end > max_coord[chrom]:
            max_coord[chrom] = end

    # output 1 record per chrom
    print_reads_intervals = [f"{chrom}\t{min_coord[chrom]}\t{max_coord[chrom]}" for chrom in max_coord]

else:
    print_reads_intervals = [f"{chrom}\t{start}\t{end}" for chrom, start, end in records]

print(f"compute_print_reads_intervals: outputing {len(print_reads_intervals)} intervals.", file=sys.stderr)
for line in print_reads_intervals:
    print(line, file=sys.stdout)

print("compute_print_reads_intervals: Done", file=sys.stderr)
