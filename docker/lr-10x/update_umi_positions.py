import subprocess

import pysam
import argparse
from Bio.Seq import Seq
import time
import os.path
import gzip
import ctypes
import sys
import shutil

from collections import defaultdict

import struct

import numpy as np

from tqdm import tqdm
from functools import reduce
import operator
import array

sys.path.append('/lrma')
from ssw import ssw_lib

NUM_BASES_TO_ALIGN_CCS = 80
NUM_BASES_TO_ALIGN_CLR = 120

UMI_LENGTH = 10
CELL_BARCODE_SEQUENCE = "TCTACACGACGCTCTTCCGATCT"
POST_UMI_SEQ = "TTTCTTATATGGG"
CELL_BARCODE_TAG = "CB"

UMI_TAG = "ZU"
NEW_UMI_TAG = "ZX"


# IUPAC RC's from: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
# and https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
RC_BASE_MAP = {
    "N": "N",
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "Y": "R",
    "R": "Y",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "n": "n",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "y": "r",
    "r": "y",
    "s": "s",
    "w": "w",
    "k": "m",
    "m": "k",
    "b": "v",
    "v": "b",
    "d": "h",
    "h": "d",
}


def reverse_complement(base_string):
    """
    Reverse complements the given base_string.
    :param base_string: String of bases to be reverse-complemented.
    :return: The reverse complement of the given base string.
    """

    return "".join(map(lambda b: RC_BASE_MAP[b], base_string[::-1]))


def ssw_build_matrix(match_score=2, mismatch_score=1):
    """
    Builds a matrix of ctypes values for matches and mismatches for the Smith-Waterman algorithm
    :param match_score: The score for matches
    :param mismatch_score: The penalty for mismatches, i.e. a positive number results in a penalty of that number.
    :return: A matrix for use in the ssw algorithm
    """
    alphabet = ['A', 'C', 'G', 'T', 'N']
    # dRc = {'A': 'C', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'C', 'c': 'G', 'g': 'C', 't': 'A'}
    letter_to_int = dict()
    int_to_letter = dict()
    for i, letter in enumerate(alphabet):
        letter_to_int[letter] = i
        letter_to_int[letter.lower()] = i
        int_to_letter[i] = letter
    alphabet_size = len(alphabet)
    lScore = [0 for i in range(alphabet_size ** 2)]
    for i in range(alphabet_size - 1):
        for j in range(alphabet_size - 1):
            if alphabet[i] == alphabet[j]:
                lScore[i * alphabet_size + j] = match_score
            else:
                lScore[i * alphabet_size + j] = -mismatch_score

    mat = (ctypes.c_int8 * len(lScore))(*lScore)
    return alphabet, letter_to_int, mat


def to_int(seq, lEle, dEle2Int):
    """
    Translates a sequence into ctypes numbers
    @param  seq   a sequence
    """
    num_decl = len(seq) * ctypes.c_int8
    num = num_decl()
    for i, ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n

    return num


def get_alignment(ssw, sequence, adapter_seq, alphabet, letter_to_int, mat, open_penalty=2, extension_penalty=1):
    """
    Performs the alignment of the read end to the adapter sequence
    :param ssw: ssw object for performing the Smith-Waterman alignment
    :param sequence: The sequence of the read end
    :param adapter_seq: The adapter sequence to align to
    :param alphabet: The alphabet object for the ssw algorithm
    :param letter_to_int: The dict to convert base letters to numbers
    :param mat: The match-mismatch score matrix for the ssw algorithm
    :param open_penalty: The penalty for opening gaps in the alignment
    :param extension_penalty: The penalty for extending gaps
    :return: A tuple containing the position of the first and last base of the adapter sequence in the read end
    """
    sequence_numbers = to_int(sequence, alphabet, letter_to_int)
    adapter_numbers = to_int(adapter_seq, alphabet, letter_to_int)

    flag = 1
    mask_length = len(adapter_seq) // 2 if len(adapter_seq) >= 30 else 15

    q_profile = ssw.ssw_init(sequence_numbers, ctypes.c_int32(len(sequence)), mat, len(alphabet), 2)

    res = ssw.ssw_align(q_profile, adapter_numbers, ctypes.c_int32(len(adapter_seq)), open_penalty,
                        extension_penalty, flag, 0, 0, mask_length)

    return res.contents.nQryBeg, res.contents.nQryEnd, res.contents.nScore, res.contents


def print_alignment(seq, adapter_seq, query_start, query_end, ref_start, ref_end, post_alignment_start, umi_length, cigar_length, cigar):

    # Get our cigar operations here:
    num_insertions = 0
    num_deletions = 0
    cigar_ops = []
    for i in range(cigar_length):
        count, op = unpack_cigar_count_and_op(cigar[i])
        if op == "I":
            num_insertions += count
        elif op == "D":
            num_deletions += count
        cigar_ops.append((count, op))

    extra_read_padding = " " * (ref_start - query_start)

    print()
    print(extra_read_padding, end="")
    print(" " * query_start, end="")
    print("V", end="")
    print("-" * (query_end - query_start - 1 + num_deletions), end="")
    print("|", end="")
    print("V", end="")
    print("-" * (umi_length - 2), end="")
    print("|")
    print(extra_read_padding, end="")
    print(" " * (post_alignment_start + num_deletions), end="")
    print("V")

    # Print the excerpt from the read:
    print(extra_read_padding, end="")
    print(seq[:query_start], end="")
    cur_pos = query_start
    for count, op in cigar_ops:
        if op == "M" or op == "I":
            print(seq[cur_pos:cur_pos+count], end="")
            cur_pos += count
        elif op == "D":
            print("-" * count, end="")
    print(seq[cur_pos:], end="")
    print()

    print(" " * (query_start - ref_start), end="")

    # Print the excerpt from the read:
    marker_len = 0
    cur_pos = 0
    for count, op in cigar_ops:
        if op == "M" or op == "D":
            print(adapter_seq[cur_pos:cur_pos+count], end="")
            cur_pos += count
        elif op == "I":
            print("-" * count, end="")
            marker_len += count
    print(adapter_seq[cur_pos:], end="")
    marker_len += cur_pos
    print(" " * ((post_alignment_start + num_deletions) - marker_len - (query_start - ref_start) - ref_start), end="")
    print(POST_UMI_SEQ)
    print()

    print(" " * (query_start - ref_start), end="")
    print(" " * ref_start, end="")
    print("^", end="")
    print("-" * (marker_len - 2), end="")
    print("|", end="")
    print("^", end="")
    print("-" * (umi_length - 2), end="")
    print("|")


def unpack_cigar_count_and_op(cigar_int):

    # top 28 bits are count:
    count = (cigar_int & 0xFFFFFFF0) >> 4

    # bottom 4 bits are type:
    t = (cigar_int & 0xF)
    if t == 0:
        op = "M"
    elif t == 1:
        op = "I"
    elif t == 2:
        op = "D"
    else:
        raise RuntimeError(f"Unpacked unknown type of cigar operation: {t} from raw cigar: {cigar_int}")

    return count, op


def get_short_reads_umis(short_reads_umi_file):
    with open(short_reads_umi_file, 'r') as f:
        umis = set()
        for line in f:
            if line.startswith("#") or line.startswith("UMI"):
                continue
            umis.add(line.strip().split()[0].strip())
    return umis


def main(bam_filename, out_file_name, barcode_seq, cell_barcode_tag, umi_length, new_umi_tag, existing_umi_tag,
         short_reads_umi_file, ssw_path):

    # Read in the short reads umis:
    short_read_umis = get_short_reads_umis(short_reads_umi_file)

    # Set up our SSW objects:
    ssw = ssw_lib.CSsw(ssw_path)
    alphabet, letter_to_int, mat = ssw_build_matrix()

    # silence message about the .bai file not being found
    pysam.set_verbosity(0)

    num_reads = 0
    num_new_umis_same_as_old = 0
    num_ccs_reads_not_same_as_old = 0
    num_clr_reads_not_same_as_old = 0

    num_old_umis_in_short_reads_umis = 0
    num_old_ccs_umis_in_short_reads_umis = 0
    num_old_clr_umis_in_short_reads_umis = 0
    num_old_umis_not_in_short_reads_umis = 0
    num_old_ccs_umis_not_in_short_reads_umis = 0
    num_old_clr_umis_not_in_short_reads_umis = 0

    num_new_umis_in_short_reads_umis = 0
    num_new_ccs_umis_in_short_reads_umis = 0
    num_new_clr_umis_in_short_reads_umis = 0
    num_new_umis_not_in_short_reads_umis = 0
    num_new_ccs_umis_not_in_short_reads_umis = 0
    num_new_clr_umis_not_in_short_reads_umis = 0

    new_tenx_start_offsets = defaultdict(int)
    new_cbc_end_offsets = defaultdict(int)
    new_post_adapter_start_offsets = defaultdict(int)

    new_umi_length_hist = defaultdict(int)
    new_ccs_umi_length_hist = defaultdict(int)
    new_clr_umi_length_hist = defaultdict(int)

    num_reads_missing_venus = 0
    num_ccs_reads_missing_venus = 0
    num_clr_reads_missing_venus = 0

    num_reads_missing_cbc = 0
    num_ccs_reads_missing_cbc = 0
    num_clr_reads_missing_cbc = 0

    num_reads_missing_boreas = 0
    num_ccs_reads_missing_boreas = 0
    num_clr_reads_missing_boreas = 0

    with pysam.AlignmentFile(bam_filename, 'rb', check_sq=False, require_index=False) as bam_file, \
            tqdm(desc=f"Processing reads", unit="read", file=sys.stderr) as pbar:
        with pysam.AlignmentFile(out_file_name, 'wb', header=bam_file.header) as out_bam_file:

            for read in bam_file:

                # Parse our segments so we can use them later for stats:
                segments = dict()
                for s in read.get_tag("SG").strip().split(','):
                    name, bounds = s.split(":")
                    st, nd = bounds.split("-")
                    st = int(st)
                    nd = int(nd)
                    segments[name] = (st, nd)

                num_reads += 1

                # We grab a different number of bases if it's CCS vs CLR because CLR reads can be VERY noisy:
                is_ccs = read.get_tag("rq") > 0
                num_read_bases_to_align = NUM_BASES_TO_ALIGN_CCS if is_ccs else NUM_BASES_TO_ALIGN_CLR

                adapter_seq = barcode_seq + read.get_tag(cell_barcode_tag)

                # Let's align the CBC and UMI in the forward direction:
                alignments = get_alignment(ssw,
                                           read.query_sequence[:num_read_bases_to_align],
                                           adapter_seq,
                                           alphabet,
                                           letter_to_int,
                                           mat)

                # Now do it in the reverse direction:
                rc_seq = reverse_complement(read.query_sequence)
                rc_alignments = get_alignment(ssw,
                                              rc_seq[:num_read_bases_to_align],
                                              adapter_seq,
                                              alphabet,
                                              letter_to_int,
                                              mat)

                # Take the direction with the best score:
                if alignments[2] > rc_alignments[2]:
                    seq = read.query_sequence
                    res = alignments[3]
                    alignment_start = alignments[0]
                    alignment_end = alignments[1]

                else:
                    seq = rc_seq
                    res = rc_alignments[3]
                    alignment_start = rc_alignments[0]
                    alignment_end = rc_alignments[1]

                # Now we align the POST_UMI_SEQ to the region following our adapters
                # so we can have a bounded region within which we can grab the UMI:
                post_alignments = get_alignment(ssw,
                                                seq[alignment_end:alignment_end + 2 * (umi_length + len(POST_UMI_SEQ))],
                                                POST_UMI_SEQ,
                                                alphabet,
                                                letter_to_int,
                                                mat)

                # Adjust positions to reflect start position in string:
                post_alignment_start = post_alignments[0] + alignment_end
                post_alignments[3].nQryBeg += alignment_end + 1
                post_alignments[3].nQryEnd += alignment_end + 1

                # Get the new UMI from the read sequence:
                # new_umi_seq = seq[alignment_end + 1:alignment_end + 1 + umi_length]
                new_umi_seq = seq[alignment_end + 1:post_alignment_start]

                print()
                print("=" * 80)
                print(f"{read.query_name}")
                print()

                # Debugging:
                print_alignment(seq[:num_read_bases_to_align], adapter_seq, res.nQryBeg, res.nQryEnd, res.nRefBeg,
                                    res.nRefEnd, post_alignment_start, umi_length, res.nCigarLen, res.sCigar)

                print()
                print(f"nQryBeg = {res.nQryBeg}")
                print(f"nQryEnd = {res.nQryEnd}")
                print(f"nRefBeg = {res.nRefBeg}")
                print(f"nRefEnd = {res.nRefEnd}")
                print(f"nRefEnd2 = {res.nRefEnd2}")
                print(f"nScore = {res.nScore}")
                print(f"nScore2 = {res.nScore2}")
                print(f"nCigarLen = {res.nCigarLen}")
                print(f"sCigar = ", end="")
                for i in range(res.nCigarLen):
                    count, op = unpack_cigar_count_and_op(res.sCigar[i])
                    print(f"{count}{op} ", end="")
                print()
                if read.has_tag(existing_umi_tag) and read.get_tag(existing_umi_tag) != new_umi_seq:
                    print(f"OLD: {read.get_tag(existing_umi_tag)}")
                    print(f"NEW: {new_umi_seq}")
                    print()
                else:
                    print(f"UMI: {new_umi_seq}")
                    print()
                print(f"Read RQ: {read.get_tag('rq'):0.5f}")
                print(f"UMI Length: {len(new_umi_seq)}")

                # Get our stats:
                if read.has_tag(existing_umi_tag):
                    old_umi = read.get_tag(existing_umi_tag)

                    # Get old stats:
                    if old_umi in short_read_umis:
                        num_old_umis_in_short_reads_umis += 1
                        if is_ccs:
                            num_old_ccs_umis_in_short_reads_umis += 1
                        else:
                            num_old_clr_umis_in_short_reads_umis += 1
                    else:
                        num_old_umis_not_in_short_reads_umis += 1
                        if is_ccs:
                            num_old_ccs_umis_not_in_short_reads_umis += 1
                        else:
                            num_old_clr_umis_not_in_short_reads_umis += 1

                    # Get new stats:
                    if old_umi == new_umi_seq:
                        num_new_umis_same_as_old += 1
                    else:
                        if is_ccs:
                            num_ccs_reads_not_same_as_old += 1
                        else:
                            num_clr_reads_not_same_as_old += 1

                # Get new stats:
                if new_umi_seq in short_read_umis:
                    num_new_umis_in_short_reads_umis += 1
                    if is_ccs:
                        num_new_ccs_umis_in_short_reads_umis += 1
                    else:
                        num_new_clr_umis_in_short_reads_umis += 1
                else:
                    num_new_umis_not_in_short_reads_umis += 1
                    if is_ccs:
                        num_new_ccs_umis_not_in_short_reads_umis += 1
                    else:
                        num_new_clr_umis_not_in_short_reads_umis += 1

                try:
                    new_tenx_start_offsets[alignment_start - segments['VENUS'][0]] += 1
                    print(f"10x start offset: {alignment_start - segments['VENUS'][0]}")
                except KeyError:
                    num_reads_missing_venus += 1
                    print("Missing 10x adapter")
                    if is_ccs:
                        num_ccs_reads_missing_venus += 1
                    else:
                        num_clr_reads_missing_venus += 1

                try:
                    new_cbc_end_offsets[alignment_end - segments['CBC'][1]] += 1
                    print(f"CBC end offset: {alignment_end - segments['CBC'][1]}")
                except KeyError:
                    num_reads_missing_cbc += 1
                    print("Missing CBC")
                    if is_ccs:
                        num_ccs_reads_missing_cbc += 1
                    else:
                        num_clr_reads_missing_cbc += 1

                try:
                    new_post_adapter_start_offsets[post_alignment_start - segments['BOREAS'][0]] += 1
                    print(f"Post UMI start offset: {post_alignment_start - segments['BOREAS'][0]}")
                except KeyError:
                    num_reads_missing_boreas += 1
                    print("Missing Post UMI Adapter")
                    if is_ccs:
                        num_ccs_reads_missing_boreas += 1
                    else:
                        num_clr_reads_missing_boreas += 1

                new_umi_length_hist[len(new_umi_seq)] += 1
                if is_ccs:
                    new_ccs_umi_length_hist[len(new_umi_seq)] += 1
                else:
                    new_clr_umi_length_hist[len(new_umi_seq)] += 1

                print()

                # Actually set the new tag for the new UMI and write the read:
                read.set_tag(new_umi_tag, new_umi_seq)

                out_bam_file.write(read)
                pbar.update(1)

    print("=" * 80)
    print()
    print("Stats:")
    print(f"Total Num Reads: {num_reads}")
    print()
    print(f"Num new UMIs same as old UMIs: {num_new_umis_same_as_old} ({100*num_new_umis_same_as_old/num_reads:2.4f}%)")
    print(f"Num new UMIs different from old UMIs: {num_reads - num_new_umis_same_as_old} ({100*(num_reads - num_new_umis_same_as_old)/num_reads:2.4f}%)")
    print(f"Num CCS Reads with new UMIs != old UMIs: {num_ccs_reads_not_same_as_old} ({100*num_ccs_reads_not_same_as_old/num_reads:2.4f}%)")
    print(f"Num CLR Reads with new UMIs != old UMIs: {num_clr_reads_not_same_as_old} ({100*num_clr_reads_not_same_as_old/num_reads:2.4f}%)")
    print()
    print(f"Num old umis in short reads umis: {num_old_umis_in_short_reads_umis} ({100*num_old_umis_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num old CCS umis in short reads umis: {num_old_ccs_umis_in_short_reads_umis} ({100*num_old_ccs_umis_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num old CLR umis in short reads umis: {num_old_clr_umis_in_short_reads_umis} ({100*num_old_clr_umis_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num old umis NOT in short reads umis: {num_old_umis_not_in_short_reads_umis} ({100*num_old_umis_not_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num old CCS umis NOT in short reads umis: {num_old_ccs_umis_not_in_short_reads_umis} ({100*num_old_ccs_umis_not_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num old CLR umis NOT in short reads umis: {num_old_clr_umis_not_in_short_reads_umis} ({100*num_old_clr_umis_not_in_short_reads_umis/num_reads:2.4f}%)")
    print()
    print(f"Num new umis in short reads umis: {num_new_umis_in_short_reads_umis} ({100*num_new_umis_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num new CCS umis in short reads umis: {num_new_ccs_umis_in_short_reads_umis} ({100*num_new_ccs_umis_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num new CLR umis in short reads umis: {num_new_clr_umis_in_short_reads_umis} ({100*num_new_clr_umis_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num new umis NOT in short reads umis: {num_new_umis_not_in_short_reads_umis} ({100*num_new_umis_not_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num new CCS umis NOT in short reads umis: {num_new_ccs_umis_not_in_short_reads_umis} ({100*num_new_ccs_umis_not_in_short_reads_umis/num_reads:2.4f}%)")
    print(f"Num new CLR umis NOT in short reads umis: {num_new_clr_umis_not_in_short_reads_umis} ({100*num_new_clr_umis_not_in_short_reads_umis/num_reads:2.4f}%)")
    print()
    print(f"Num reads missing 10x adapter (venus): {num_reads_missing_venus} ({100 * num_reads_missing_venus / num_reads:2.4f}%)")
    print(f"Num CCS reads missing 10x adapter (venus): {num_ccs_reads_missing_venus} ({100 * num_ccs_reads_missing_venus / num_reads:2.4f}%)")
    print(f"Num CLR reads missing 10x adapter (venus): {num_clr_reads_missing_venus} ({100 * num_clr_reads_missing_venus / num_reads:2.4f}%)")
    print()
    print(f"Num reads missing CBC: {num_reads_missing_cbc} ({100 * num_reads_missing_cbc / num_reads:2.4f}%)")
    print(f"Num CCS reads missing CBC: {num_ccs_reads_missing_cbc} ({100 * num_ccs_reads_missing_cbc / num_reads:2.4f}%)")
    print(f"Num CLR reads missing CBC: {num_clr_reads_missing_cbc} ({100 * num_clr_reads_missing_cbc / num_reads:2.4f}%)")
    print()
    print(f"Num reads missing post-umi adapter (boreas): {num_reads_missing_boreas} ({100 * num_reads_missing_boreas / num_reads:2.4f}%)")
    print(f"Num CCS reads missing post-umi adapter (boreas): {num_ccs_reads_missing_boreas} ({100 * num_ccs_reads_missing_boreas / num_reads:2.4f}%)")
    print(f"Num CLR reads missing post-umi adapter (boreas): {num_clr_reads_missing_boreas} ({100 * num_clr_reads_missing_boreas / num_reads:2.4f}%)")
    print()
    print("Overall new UMI Lengths:")
    print_dict_stats(new_umi_length_hist)
    print()
    print("CCS new UMI Lengths:")
    print_dict_stats(new_ccs_umi_length_hist)
    print()
    print("CLR new UMI Lengths:")
    print_dict_stats(new_clr_umi_length_hist)
    print()
    print("New 10x adapter start offset stats:")
    print_dict_stats(new_tenx_start_offsets)

    print("New CBC end offset stats:")
    print_dict_stats(new_cbc_end_offsets)

    print("New post UMI adapter start offset stats:")
    print_dict_stats(new_post_adapter_start_offsets)
    print()
    print()
    print("For Offset histograms:")
    print("  A negative number indicates that the new offset is before the original offset.")
    print("  A positive number indicates that the new offset is after the original offset.")
    print()


def print_dict_stats(d, front_padding=2):

    d = {k:d[k] for k in sorted(d.keys())}

    print((front_padding * ' ') + "Histogram:")
    for k, v in d.items():
        print((front_padding * ' ') + f"  {k}  {v}")

    # reconstruct the base data array for ease of statistical calculations:
    data = []
    for k, v in d.items():
        for i in range(v):
            data.append(k)
    data = np.array(data)

    print((front_padding * ' ') + f"Min {np.min(data)}")
    print((front_padding * ' ') + f"Max {np.max(data)}")
    print((front_padding * ' ') + f"Mean {np.mean(data)}")
    print((front_padding * ' ') + f"Median {np.median(data)}")
    print((front_padding * ' ') + f"Stdev {np.std(data)}")
    print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Align 10x adapter + cell barcode to front of read, get next 10 bases as putative barcode.",
        epilog="")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument('-b', '--bam', help='BAM filename', required=True)
    requiredNamed.add_argument('-s', '--short-reads-umis', help='File containing short reads UMIs', required=True)
    requiredNamed.add_argument('-o', "--output-file", help='Output file', required=True)

    parser.add_argument('--ssw-path', help='Path to the Striped Smith-Waterman library', type=str, default='/lrma/ssw')

    # User-defined length for the marker segments:
    parser.add_argument(
        '--barcode-seq', help='Sequence of the cell barcode used in the library prep.  (default: %(default)s)', type=str,
        default=CELL_BARCODE_SEQUENCE
    )
    parser.add_argument(
        '--cell-barcode-tag', help='Read tag containing the corrected cell barcode for each read.  (default: %(default)s)', type=str,
        default=CELL_BARCODE_TAG
    )
    parser.add_argument(
        '--umi-length', help='Length of the unique molecular identifier used in the library prep.  (default: %(default)d)', type=int,
        default=UMI_LENGTH
    )
    parser.add_argument(
        '--new-umi-tag', help='Read tag to populate with the UMIs discovered by this tool.  (default: %(default)s)', type=str,
        default=NEW_UMI_TAG
    )
    parser.add_argument(
        '--existing-umi-tag', help='Read tag of the existing UMI in the reads.  (default: %(default)s)', type=str,
        default=UMI_TAG
    )

    args = parser.parse_args()

    main(args.bam, args.output_file, args.barcode_seq, args.cell_barcode_tag, args.umi_length, args.new_umi_tag,
         args.existing_umi_tag, args.short_reads_umis, args.ssw_path)
