#!/usr/bin/env python

import sys
import os
import re
import pysam
from tqdm import tqdm

if len(sys.argv) != 3:
    print(f"ERROR: You must give exactly 2 arguments: BAM, PREFIX", file=sys.stderr)
    sys.exit(1)

bam = sys.argv[1]
prefix = sys.argv[2]

PADDING = 2
SEG_TAG = "SG"
SEG_DELIM = ","

pysam.set_verbosity(0)

def get_read_count_from_bam_index(bam_file_path):
    """Return the number of reads in the given bam file if a bam index is present.

    Input should be a file-like object. If no index is present or input is not seekable, returns `None`.
    """

    with pysam.AlignmentFile(bam_file_path, "rb", check_sq=False, require_index=False) as bam_file:
        # Get total number of reads if we have an index:
        if bam_file.has_index():
            idx_stats = bam_file.get_index_statistics()
            unaligned_reads = bam_file.nocoordinate
            aligned_reads = reduce(lambda a, b: a + b, [x.total for x in idx_stats]) if len(idx_stats) > 0 else 0
            total_reads = unaligned_reads + aligned_reads

    return total_reads

def load_read_count(pbi_file):
    """Compute file offsets for specified read names"""

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need
    # until we get to the read count.
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
        )

    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        return idx_contents.n_reads

################################################################################

pbi = f"{bam}.pbi"
num_reads = load_read_count(pbi) if os.path.exists(pbi) else None
if not num_reads:
    num_reads = get_read_count_from_bam_index(bam)
if num_reads:
    print(f"Tot reads: {num_reads}")


with pysam.AlignmentFile(f"{bam}", "rb", check_sq=False, require_index=False) as bam_file:

    with pysam.AlignmentFile(f"{prefix}.bam", "wb", header=bam_file.header) as out_bam_file:
        for read in tqdm(bam_file, total=num_reads):
            # Fix seg coords:
            if read.is_reverse:
                segments_tag = read.get_tag(SEG_TAG)
                read_len = len(read.query_sequence)
                segments = segments_tag.split(SEG_DELIM)
                new_segments = []
                for seg in segments:
                    seg_name, start, end = re.split("[:-]", seg)
                    new_start = read_length - int(end)
                    new_end = read_length - int(start)
                    new_segments.append(f"{seg_name}:{new_start}-{new_end}")
                read.set_tag(SEG_TAG, SEG_DELIM.join(new_segments))

            # Fix original name:
            for s in read.get_tag(SEG_TAG).split(SEG_DELIM):
                seg_name, start, end = re.split("[:-]", s)
                if seg_name == "cDNA":
                    read.query_name = read.query_name + f"/{int(start)-PADDING}_{int(end)+PADDING}"
                    break

            out_bam_file.write(read)


