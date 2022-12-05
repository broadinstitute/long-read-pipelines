#!/usr/bin/env python

import sys
import pysam
import ctypes
import pickle
import argparse

import editdistance
import numpy as np
from tqdm import tqdm
from functools import reduce

from ssw import ssw_lib
ssw_path = "/usr/local/lib/ssw"

###################################

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
    dRc = {'A': 'C', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'C', 'c': 'G', 'g': 'C', 't': 'A'}
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

def get_alignment(ssw, target_seq, ref_seq, alphabet, letter_to_int, mat, open_penalty=2, extension_penalty=1):
    """
    Performs the alignment of the read end to the adapter sequence
    :param ssw: ssw object for performing the Smith-Waterman alignment
    :param sequence: The sequence of the read end
    :param adapter_sequence: The adapter sequence to align to
    :param alphabet: The alphabet object for the ssw algorithm
    :param letter_to_int: The dict to convert base letters to numbers
    :param mat: The match-mismatch score matrix for the ssw algorithm
    :param open_penalty: The penalty for opening gaps in the alignment
    :param extension_penalty: The penalty for extending gaps
    :return: A tuple containing the position of the first and last base of the adapter sequence in the read end
    """
    target_numbers = to_int(target_seq, alphabet, letter_to_int)
    ref_numbers = to_int(ref_seq, alphabet, letter_to_int)

    flag = 1
    mask_length = len(ref_seq) // 2 if len(ref_seq) >= 30 else 15

    q_profile = ssw.ssw_init(target_numbers, ctypes.c_int32(len(target_seq)), mat, len(alphabet), 2)

    res = ssw.ssw_align(q_profile, ref_numbers, ctypes.c_int32(len(ref_seq)), open_penalty,
                        extension_penalty, flag, 0, 0, mask_length)

    return res

###################################

ssw = ssw_lib.CSsw(ssw_path)

def get_cigar_tuples(cigar, cigar_len):
    cigar_tuples = []
    for i in range(cigar_len):
        op_len = (cigar[i] & 0xFFFFFFF0) >> 4
        op_char = cigar[i] & 0x0000000F

        if op_char == 0:
            op_char = "M"
        elif op_char == 1:
            op_char = "I"
        else:
            op_char = "D"
        
        cigar_tuples.append((op_len, op_char))
    return cigar_tuples

alphabet, letter_to_int, mat = ssw_build_matrix()

ssw_align = lambda x, y: get_alignment(ssw, x, y, alphabet, letter_to_int, mat)

################################################################################

# We actually need an arg parser now:
parser = argparse.ArgumentParser(
    description=f"Extracts the middle 10 bases of the 3' adapter to a tag which can then later be used for a UMI analysis.",
)

required_named_args = parser.add_argument_group('required named arguments')
required_named_args.add_argument('-b', '--bam',
                                 help='Annotated bam file from which to create the counts matrix.',
                                 required=True)
required_named_args.add_argument('-o', '--out-prefix',
                                 help='Prefix of the output files.',
                                 required=True)

optional_named_args = parser.add_argument_group("optional named arguments")
optional_named_args.add_argument('-s', '--min-sw-score',
                                 type=int,
                                 help=f"Minimum Smith-Waterman alignment score for the 3' region in the read to be "
                                      f"considered good enough for inclusion in the analysis.  Reads below this score are "
                                      f"output in a separate file.",
                                 default=35,
                                 required=False)

# Parse our args:
args = parser.parse_args()
bam = args.bam
prefix = args.out_prefix
FULL_ADAPTER_SSW_SCORE_THRESH = args.min_sw_score

# Constants:
THREE_P_ADAPTER_REAL_SEQ = "GTACTCTGCGTTGATACCACTGCTT"
THREE_P_ADAPTER_MIDDLE = THREE_P_ADAPTER_REAL_SEQ[7:17]  # GCGTTGATAC

TEST_BARCODE_TAG = "ZV"

three_p_adapter_tag = "Ja"
three_p_adapter_lev_tag = "Jb"
three_p_adapter_sw_tag = "Jc"
three_p_adapter_sw2_tag = "Jd"

num_no_three_p_adapter = 0
num_many_three_p_adapter = 0

num_no_threep = 0

padding = 15

ccs_levs = []
clr_levs = []

ccs_sws = []  
clr_sws = []

# Do the work:
with pysam.AlignmentFile(f"{bam}", "rb", check_sq=False, require_index=False) as bam_file:

    total_reads = None 
    if bam_file.has_index():
        idx_stats = bam_file.get_index_statistics()
        unaligned_reads = bam_file.nocoordinate
        aligned_reads = reduce(lambda a, b: a + b, [x.total for x in idx_stats]) if len(idx_stats) > 0 else 0
        total_reads = unaligned_reads + aligned_reads

    with pysam.AlignmentFile(f"{prefix}.bam", "wb", header=bam_file.header) as out_bam_file:
        with pysam.AlignmentFile(f"{prefix}.rejected_ssw_score_below_{FULL_ADAPTER_SSW_SCORE_THRESH}.bam", "wb", header=bam_file.header) as rejected_bam_file:
            with pysam.AlignmentFile(f"{prefix}.rejected_no_threep.bam", "wb", header=bam_file.header) as rejected_no_threep_bam_file:
                for read in tqdm(bam_file, desc="Processing reads", total=total_reads, unit=" read"):

                    three_p_adapter_seg = [s for s in read.get_tag("SG").split(",") if s.startswith("3p_Adapter")]
                    if len(three_p_adapter_seg) == 1:
                        # Get the coordinates:
                        start, end = three_p_adapter_seg[0].split(":")[1].split("-")
                        start = int(start) - padding
                        end = int(end) + padding

                        if start < 0:
                            start = 0
                        if end > len(read.query_sequence) - 1:
                            end = len(read.query_sequence) - 1

                        if read.is_reverse:
                            three_p_adapter_seq = reverse_complement(read.query_sequence)[start:end+1]
                        else:
                            three_p_adapter_seq = read.query_sequence[start:end+1]

                        res = ssw_align(three_p_adapter_seq, THREE_P_ADAPTER_REAL_SEQ)
                        read.set_tag(f"{three_p_adapter_sw_tag}", res.contents.nScore)
                        read.set_tag(f"{three_p_adapter_sw2_tag}", res.contents.nScore2)
                        three_p_adapter_in_read = three_p_adapter_seq[res.contents.nQryBeg:res.contents.nQryEnd+1]
                        
                        read.set_tag(f"{three_p_adapter_tag}", three_p_adapter_in_read)

                        lev = editdistance.eval(THREE_P_ADAPTER_REAL_SEQ, three_p_adapter_in_read)
                        read.set_tag(f"{three_p_adapter_lev_tag}", lev)

                        full_adapter_alignment_score = res.contents.nScore

                        # Now get the alignment to the MIDDLE of the 3p adapter so we can use it as a UMI:
                        res = ssw_align(three_p_adapter_in_read, THREE_P_ADAPTER_MIDDLE)
                        threep_extract = three_p_adapter_in_read[res.contents.nQryBeg:res.contents.nQryEnd+1]

                        if full_adapter_alignment_score >= FULL_ADAPTER_SSW_SCORE_THRESH:
                            if threep_extract:
                                read.set_tag(f"{TEST_BARCODE_TAG}", threep_extract)
                                out_bam_file.write(read)
                            else:
                                rejected_no_threep_bam_file.write(read)
                                num_no_threep += 0
                        else:
                            rejected_bam_file.write(read)

                        if read.get_tag("rq") >= 0.99:
                            ccs_levs.append(lev)
                            ccs_sws.append(res.contents.nScore)
                        else:
                            clr_levs.append(lev)
                            clr_sws.append(res.contents.nScore)

                    elif len(three_p_adapter_seg) == 0:
                        num_no_three_p_adapter += 1
                    else:
                        num_many_three_p_adapter += 1

ccs_levs = np.array(ccs_levs)
clr_levs = np.array(clr_levs)

ccs_sws = np.array(ccs_sws)
clr_sws = np.array(clr_sws)

pickle.dump(ccs_levs, open(f"{prefix}.ccs_levs.pickle", 'wb'))
pickle.dump(clr_levs, open(f"{prefix}.clr_levs.pickle", 'wb')) 

pickle.dump(ccs_sws, open(f"{prefix}.ccs_sws.pickle", 'wb'))
pickle.dump(clr_sws, open(f"{prefix}.clr_sws.pickle", 'wb')) 

bin_width = 1
bin_edges = np.arange(-.5, 50, bin_width)
plotty_bins = bin_edges[:-1] + bin_width/2

ccs_hist, _ = np.histogram(ccs_levs, bin_edges)
clr_hist, _ = np.histogram(clr_levs, bin_edges)

# Print table:
print("CCS")
print(f"Lev Dist     Count")
print("-"*40)
num_places = int(np.log10(np.max(ccs_hist)))+1
for l, c in zip(plotty_bins, ccs_hist):
    if c > 0:
        print(f"{int(l):8d}     {c:{num_places}d}")
        
print("CLR")
print(f"Lev Dist     Count")
print("-"*40)
num_places = int(np.log10(np.max(clr_hist)))+1
for l, c in zip(plotty_bins, clr_hist):
    if c > 0:
        print(f"{int(l):8d}     {c:{num_places}d}")

################################################################################

print()
print()

bin_width = 1
bin_edges = np.arange(-.5, 51, bin_width)
plotty_bins = bin_edges[:-1] + bin_width/2

ccs_hist, _ = np.histogram(ccs_sws, bin_edges)
clr_hist, _ = np.histogram(clr_sws, bin_edges)

# Print table:
print("CCS")
print(f"SW Score     Count")
print("-"*40)
num_places = int(np.log10(np.max(ccs_hist)))+1
for l, c in zip(plotty_bins, ccs_hist):
    if c > 0:
        print(f"{int(l):8d}     {c:{num_places}d}")
        
print("CLR")
print(f"SW Score    Count")
print("-"*40)
num_places = int(np.log10(np.max(clr_hist)))+1
for l, c in zip(plotty_bins, clr_hist):
    if c > 0:
        print(f"{int(l):8d}     {c:{num_places}d}")

################################################################################

print(f"Num no three_p_adapter: {num_no_three_p_adapter}")
print(f"Num many three_p_adapter: {num_many_three_p_adapter}")

print(f"Num no threep: {num_no_threep}")
