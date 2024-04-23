#!/usr/bin/env python3

import numpy as np
import pysam
import argparse
from tqdm import tqdm

def n50(lengths):
    all_len = sorted(lengths, reverse=True)
    csum = np.cumsum(all_len)
    n2 = int(sum(lengths) / 2)
    csumn2 = min(csum[csum >= n2])
    ind = np.where(csum == csumn2)[0]

    return all_len[int(ind[0])]


def get_bam_stats(bam_file_path, qual_thresh=None):
    # Open the file and get ready to iterate:
    with pysam.AlignmentFile(bam_file_path, "rb", check_sq=False, require_index=False) as bam_file:

        # Get total number of reads if we have an index:
        total_reads = None
        if bam_file.has_index():
            idx_stats = bam_file.get_index_statistics()
            unaligned_reads = bam_file.nocoordinate
            aligned_reads = reduce(lambda a, b: a + b, [x.total for x in idx_stats]) if len(idx_stats) > 0 else 0
            total_reads = unaligned_reads + aligned_reads

        n_reads = 0 if not total_reads else total_reads
        read_lengths = []
        quals = []
        total_bases = 0

        # Iterate through our reads
        for read in tqdm(bam_file, desc=f"Collecting Bam Stats" + (f" (rq >= {qual_thresh})" if qual_thresh else ""),
                         total=total_reads, unit=" read"):
            l = len(read.query_sequence)
            if read.query_qualities is not None:
                q = np.mean(read.query_qualities)
            else:
                q = 0

            if qual_thresh and q < qual_thresh:
                continue

            quals.append(q)
            total_bases += l
            read_lengths.append(l)

            if not total_reads:
                n_reads += 1

    return n_reads, total_bases, np.mean(quals), np.median(quals), np.array(read_lengths)


def nan_zero_wrap(float_val: float) -> float:
    """Return 0 if the value is NaN, otherwise return the value."""
    return float_val if not np.isnan(float_val) else 0


def main():
    parser = argparse.ArgumentParser(description='Compute short read bam file stats', prog='compute_sr_stats')
    parser.add_argument('-q', '--qual-threshold', type=int, default=0, help="Phred-scale quality threshold")
    parser.add_argument('bam_file_path', type=str, help="Path to bam file")
    args = parser.parse_args()

    n_reads, n_bases, mean_qual, median_qual, read_lengths = get_bam_stats(args.bam_file_path, args.qual_threshold)

    print(f"reads\t{nan_zero_wrap(n_reads)}")
    print(f"bases\t{nan_zero_wrap(n_bases)}")
    print(f"mean_qual\t{nan_zero_wrap(mean_qual)}")
    print(f"median_qual\t{nan_zero_wrap(median_qual)}")

    print(f"read_mean\t{int(np.mean(read_lengths)) if len(read_lengths) > 0 else 0}")
    print(f"read_median\t{int(np.median(read_lengths)) if len(read_lengths) > 0 else 0}")
    print(f"read_stdev\t{int(np.std(read_lengths)) if len(read_lengths) > 0 else 0}")
    print(f"read_n50\t{n50(read_lengths) if len(read_lengths) > 0 else 0}")


if __name__ == "__main__":
    main()

