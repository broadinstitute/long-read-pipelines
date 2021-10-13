#!/usr/bin/env python

import pysam
import time
import math
import sys
import os
import ssw
import gzip

from construct import *

from tqdm import tqdm

index_sequences = [
    "TTCTTGGCGG",  # 1
    "AGCGGATCGA",  # 2
    "TGTCACCATT",  # 3
    "TTCTCCTGAA",  # 4
]

# base_folder = "/home/jovyan/work/juffowup2/MAS_seq/SIRV_MOD"
# input_bam = "SIRV_MASx15_mod_preps_all_annotated_reads.1000.bam"
# input_bam = "SIRV_MASx15_mod_preps_all_annotated_reads.10000.bam"

window_size = 5
ssw_aligner = ssw.Aligner()

max_iter = math.inf
cur_iter = 0

ambiguity_thresh_pct = 20
num_ambiguous = 0
num_reads_demuxed_by_index = [0 for s in index_sequences]

index_tag = "XI"

################################################################################


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

t_start = time.time()

base_folder = "."
try:
    input_bam = sys.argv[1]
except IndexError:
    print("ERROR: You must supply an input file.", file=sys.stderr)
    sys.exit(1)

################################################################################

with pysam.AlignmentFile(f"{input_bam}", "rb", check_sq=False, require_index=False) as bam_file, tqdm(
        desc="Demultiplexing indexed reads", unit="read") as pbar:

    out_base_file_name = input_bam[:input_bam.find(".bam")] + ".demux_"

    pbi = f"{input_bam}.pbi"
    read_count = None
    if os.path.exists(pbi):
        read_count = load_read_count(pbi)

    epsilon = 1 - (ambiguity_thresh_pct / 100)

    bam_header_dict = bam_file.header.to_dict()
    pg_dict = {
        "ID": f"demux_by_index",
        "PN": "demux_by_index",
        "VN": f"0.0.1",
        "DS": "Demultiplex the given bam file using the index base sequence.",
        "CL": f"demux_by_index --window-size {window_size} --ambiguity-thresh-pct {ambiguity_thresh_pct}",
    }
    if "PG" in bam_header_dict:
        bam_header_dict["PG"].append(pg_dict)
    else:
        bam_header_dict["PG"] = [pg_dict]
    out_bam_header = pysam.AlignmentHeader.from_dict(bam_header_dict)

    try:
        out_files = []
        for i in range(len(index_sequences)):
            f_name = f"{base_folder}/{out_base_file_name}i{i + 1}.bam"
            print(f"Creating output file for index {i + 1}: {f_name}")
            of = pysam.AlignmentFile(f_name, "wb", header=out_bam_header)
            out_files.append(of)
        f_name = f"{base_folder}/{out_base_file_name}ambiguous_indices.bam"
        print(f"Creating output file for ambiguous indices: {f_name}")
        out_files.append(pysam.AlignmentFile(f_name, "wb", header=out_bam_header))

        with tqdm(desc="Demultiplexing indexed reads", unit="read", total=read_count) as pbar:
            for read in bam_file:

                raw_scores = [0 for s in index_sequences]

                segments = read.get_tag("SG").split(",")
                index_containing_adapters = [s for s in segments if s.startswith('3p_Adapter')]

                for a in index_containing_adapters:
                    # Get indices:
                    i1, i2 = a.split(":")[1].split("-")

                    # Adjust for indels:
                    i1 = int(i1)
                    i2 = i1 + 10 + window_size
                    i1 -= window_size

                    # Bounds check:
                    if i1 < 0:
                        i1 = 0
                    if i2 > len(read.query_sequence) - 1:
                        i2 = len(read.query_sequence)

                    # Extract sequence:
                    raw_index_seq = read.query_sequence[i1:i2 + 1]

                    for i, index_seq in enumerate(index_sequences):
                        alignment = ssw_aligner.align(raw_index_seq.upper(), index_seq)
                        optimal_score = alignment.score
                        max_score = len(index_seq) * ssw_aligner.matrix.get_match()

                        raw_scores[i] += optimal_score

                # Get the best score:
                max_score = 0
                max_index = 0
                approximately_equal_score_pairs = []
                for i in range(len(index_sequences)):
                    if raw_scores[i] > max_score:
                        max_score = raw_scores[i]
                        max_index = i

                # Check for possible ambiguities:
                approximately_equal_score_pairs = []
                approximately_equal_score_ratios = []
                if max_score == 0:
                    for i in range(len(index_sequences)):
                        for j in range(i + 1, len(index_sequences)):
                            approximately_equal_score_pairs.append((i, j))
                            approximately_equal_score_ratios.append(0)
                else:
                    for i in range(len(index_sequences)):
                        if i == max_index:
                            continue
                        score_ratio = raw_scores[i] / max_score
                        if score_ratio > epsilon:
                            approximately_equal_score_pairs.append((i, max_index))
                            approximately_equal_score_ratios.append(score_ratio)

                # Log and output reads to the right file:
                if len(approximately_equal_score_pairs) != 0:
                    print(
                        f"Ambiguous read: {read.query_name}: [MI={max_index}] {raw_scores} - {approximately_equal_score_pairs}: {approximately_equal_score_ratios}")
                    read.set_tag(index_tag, "AMBIGUOUS")
                    out_files[len(out_files) - 1].write(read)
                    num_ambiguous += 1
                else:
                    read.set_tag(index_tag, index_sequences[max_index])
                    out_files[max_index].write(read)
                    num_reads_demuxed_by_index[max_index] += 1

                pbar.update(1)
                cur_iter += 1

                if cur_iter > max_iter:
                    break
    finally:
        for f in out_files:
            f.close()

t_end = time.time()
print(f"Total reads seen: {cur_iter}")
print(f"Num ambiguous reads: {num_ambiguous}")
print(f"Ambiguity percentage: {100 * num_ambiguous / cur_iter:.3f}%")
print(f"Reads demuxed by index: ")
print("\t".join([str(i + 1) for i in range(len(index_sequences))]))
print("\t".join([str(n) for n in num_reads_demuxed_by_index]))
print()
print(f"Time elapsed: {t_end - t_start:.2f}s")
print(f"Time per read: {(t_end - t_start) / cur_iter:.4f}s")
