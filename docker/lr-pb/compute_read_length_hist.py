from __future__ import print_function
import sys
import argparse
import gzip
import math
import numpy
from construct import *


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def to_phred_score(p):
    return 40 if p >= 1.0 else int(-10.0 * math.log10(1.0 - p))


def n50(lengths):
    all_len = sorted(lengths, reverse=True)
    csum = numpy.cumsum(all_len)
    n2 = int(sum(lengths)/2)
    csumn2 = min(csum[csum >= n2])
    ind = numpy.where(csum == csumn2)

    return all_len[int(ind[0])]


def load_index(pbi_file, qual_threshold):
    """
    Load .pbi index data
    """

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need for sharding.
    # More on index format at https://pacbiofileformats.readthedocs.io/en/9.0/PacBioBamIndex.html .

    fmt = Struct(
        # Header
        "magic" / Const(b"PBI\x01"),
        "version_patch" / Int8ul,
        "version_minor" / Int8ul,
        "version_major" / Int8ul,
        "version_empty" / Int8ul,
        "pbi_flags" / Int16ul,
        "n_reads" / Int32ul,
        "reserved" / Padding(18),

        # Basic information section (columnar format)
        "rgId" / Padding(this.n_reads * 4),
        "qStart" / Array(this.n_reads, Int32sl),
        "qEnd" / Array(this.n_reads, Int32sl),
        "holeNumber" / Array(this.n_reads, Int32sl),
        "readQual" / Array(this.n_reads, Float32l),
    )

    n_reads = 0
    polymerase_read_lengths = {}
    subread_lengths = []
    quals = []
    total_bases = 0
    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        for j in range(0, idx_contents.n_reads):
            if to_phred_score(idx_contents.readQual[j]) >= qual_threshold:
                length = idx_contents.qEnd[j] - idx_contents.qStart[j]

                # Save the polymerase and subread lengths
                if idx_contents.holeNumber[j] not in polymerase_read_lengths:
                    polymerase_read_lengths[idx_contents.holeNumber[j]] = 0

                n_reads += 1
                polymerase_read_lengths[idx_contents.holeNumber[j]] += length
                subread_lengths.append(length)
                quals.append(to_phred_score(idx_contents.readQual[j]))
                total_bases += length

    return n_reads, total_bases, numpy.mean(quals), numpy.median(quals), polymerase_read_lengths, subread_lengths


def main():
    parser = argparse.ArgumentParser(description='Compute read length distribution from .pbi', prog='compute_read_length_hist')
    parser.add_argument('-q', '--qual-threshold', type=int, default=0, help="Phred-scale quality threshold")
    parser.add_argument('pbi', type=str, help=".pbi index")
    args = parser.parse_args()

    # Decode PacBio .pbi file and determine the polymerase and subread lengths
    eprint(f"Reading index ({args.pbi}). This may take a few minutes...", flush=True)
    n_reads, n_bases, mean_qual, median_qual, polymerase_read_lengths, subread_lengths = load_index(args.pbi, args.qual_threshold)

    for subread_length in subread_lengths:
        print(f'{subread_length}')


if __name__ == "__main__":
    main()
