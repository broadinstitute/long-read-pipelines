#!/usr/bin/env python3

import argparse
import collections
import os
import re

p = argparse.ArgumentParser(description="Takes a .bed file path and splits it into the specified number of subsets so that "
    "the total genomic region covered by the intervals in each shard is roughly equal.")
p.add_argument("--num-shards", help="Create this many output .bed files. Note: this is treated as an approximate limit. "
    "Slightly more shards may be created to that each shard contains no more than a max number of intervals (= total_intervals / num_shards)"
    "and also spans no more than the max number of bases (= total_bases_spanned / num_shards).", type=int, required=True)
p.add_argument("--output-dir-suffix", default="_sharded", help="Output .bed files will be written to an output "
    " directory that's named like the input .bed file name (without '.bed') + this suffix.")
p.add_argument("bed_file_path", nargs="+", help="One or more input .bed files to split. Each will be split separately.")
args = p.parse_args()

for bed_file_path in args.bed_file_path:
    # compute num_records in bed_file_path
    total_bases_covered = 0
    total_intervals_counter = 0

    records = collections.defaultdict(list)

    with open(bed_file_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            records[chrom].append((start, end, line))
            total_intervals_counter += 1

    # make sure intervals are sorted, otherwise too many shards could be generated
    for chrom, intervals in records.items():
        if [(i[0], i[1]) for i in intervals] != [(i[0], i[1]) for i in sorted(intervals, key=lambda i: i[0])]:
            raise ValueError(f"{bed_file_path} is not sorted. It can be sorted with bedtools sort -i {bed_file_path}")

    total_bases_covered = 0
    for chrom, intervals in records.items():
        min_start = min(i[0] for i in intervals)
        max_end = max(i[1] for i in intervals)

        total_bases_covered += (max_end - min_start)

    bases_per_shard = total_bases_covered // args.num_shards
    max_intervals_per_shard = total_intervals_counter // args.num_shards

    print(f"Parsed {sum(len(intervals) for intervals in records.values())} intervals with {total_bases_covered:,} total bases spanned. "
        f"Will span {bases_per_shard:,}bp per shard, and contain no more than {max_intervals_per_shard} intervals per shard.")

    # create output dir and delete previous files
    basename = os.path.basename(bed_file_path)
    output_filename_prefix = re.sub(".bed$", "", basename)
    output_dir = f"{output_filename_prefix}{args.output_dir_suffix}"
    os.system(f"mkdir -p {output_dir}")
    os.system(f"rm {output_dir}/{output_filename_prefix}.*.bed")

    shard_i = 0
    output_file = None
    bases_covered = 0
    line_counter = 0
    intervals_per_shard_counter = 0
    min_start_in_chrom = 0
    max_end_in_chrom = 0

    for chrom, intervals in records.items():
        for start, end, line in intervals:
            if not min_start_in_chrom or start < min_start_in_chrom:
                min_start_in_chrom = start
            if not max_end_in_chrom or end > max_end_in_chrom:
                max_end_in_chrom = end

            if output_file is None or bases_covered + (max_end_in_chrom - min_start_in_chrom) > bases_per_shard or intervals_per_shard_counter > max_intervals_per_shard:
                if output_file is not None:
                    print(f"Creating shard {shard_i}: {output_file.name} with {intervals_per_shard_counter} intervals, "
                        f"spanning {bases_covered + (max_end_in_chrom - min_start_in_chrom):,}bp.")
                    output_file.close()

                output_file = open(f"{output_dir}/{output_filename_prefix}.{shard_i:02d}.bed", "w")
                bases_covered = 0
                intervals_per_shard_counter = 0
                min_start_in_chrom = 0
                max_end_in_chrom = 0
                shard_i += 1

                if shard_i > (2 * args.num_shards):
                    raise ValueError(f"Something went wrong with sharding. Too many shards created: {shard_i}.")

            intervals_per_shard_counter += 1
            line_counter += 1
            output_file.write(line)

        bases_covered += max_end_in_chrom - min_start_in_chrom
        max_end_in_chrom = 0
        min_start_in_chrom = 0



    if output_file is not None:
        print(f"Creating shard {shard_i}: {output_file.name} with {intervals_per_shard_counter} intervals, "
            f"spanning {bases_covered + (max_end_in_chrom - min_start_in_chrom):,}bp.")
        output_file.close()

    print(f"Done. Wrote {line_counter} lines total to {shard_i} shards.")
