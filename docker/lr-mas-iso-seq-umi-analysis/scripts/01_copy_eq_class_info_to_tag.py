#!/usr/bin/env python

import sys
import os
import pysam
from tqdm import tqdm

if len(sys.argv) != 6:
    print(f"ERROR: You must give exactly 5 arguments: BAM, EQ_CLASS_FILE, EQ_CLASS_TAG, GENE_TAG, PREFIX", file=sys.stderr)
    sys.exit(1)

bam = sys.argv[1]
eq_class_file = sys.argv[2]
eq_class_tag = sys.argv[3]
gene_tag = sys.argv[4]
prefix = sys.argv[5]

################################################################################

def _blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def get_num_lines_in_text_file(file_name):
    with open(file_name, "r", encoding="utf-8",errors='ignore') as f:
        return sum(bl.count("\n") for bl in _blocks(f))

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

################################################################################

print("Counting num lines in EQ class file...")
num_lines = get_num_lines_in_text_file(f"{eq_class_file}")
print(f"Num lines: {num_lines}")
with open(f"{eq_class_file}", 'r') as f:
    read_gene_dict = dict()
    for line in tqdm(f, desc="Reading in EQ classes", total=num_lines):
        if line.startswith("#"):
            continue
        read_name, tx_eq_class, gene_assignment = line.strip().split("\t")
        read_gene_dict[read_name] = (tx_eq_class, gene_assignment)

################################################################################

with pysam.AlignmentFile(f"{bam}", "rb", check_sq=False, require_index=False) as bam_file:
    with pysam.AlignmentFile(f"{prefix}.bam", "wb", header=bam_file.header) as out_bam_file:
        for read in tqdm(bam_file, desc="Assigning EQ classes to reads", total=num_reads):
            # Not all reads from the raw data make it into the final output, so we expect
            # some reads are going to be missing:
            try:
                read.set_tag(f"{eq_class_tag}", read_gene_dict[read.query_name][0])
                read.set_tag(f"{gene_tag}", read_gene_dict[read.query_name][1])
                out_bam_file.write(read)
            except KeyError:
                pass


