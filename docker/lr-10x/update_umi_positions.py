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

import numpy as np

from tqdm import tqdm
from functools import reduce
import operator
import array

sys.path.append('/lrma')
from ssw import ssw_lib

NUM_BASES_TO_ALIGN = 80

UMI_LENGTH = 10
CELL_BARCODE_SEQUENCE = "TCTACACGACGCTCTTCCGATCT"
CELL_BARCODE_TAG = "CB"

UMI_TAG = "ZU"
NEW_UMI_TAG = "ZX"


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


def get_alignment(ssw, sequence, adapter_sequence, alphabet, letter_to_int, mat, open_penalty=2, extension_penalty=1):
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
    sequence_numbers = to_int(sequence, alphabet, letter_to_int)
    adapter_numbers = to_int(adapter_sequence, alphabet, letter_to_int)

    flag = 1
    mask_length = len(adapter_sequence) // 2 if len(adapter_sequence) >= 30 else 15

    q_profile = ssw.ssw_init(sequence_numbers, ctypes.c_int32(len(sequence)), mat, len(alphabet), 2)

    res = ssw.ssw_align(q_profile, adapter_numbers, ctypes.c_int32(len(adapter_sequence)), open_penalty,
                        extension_penalty, flag, 0, 0, mask_length)

    return (res.contents.nQryBeg, res.contents.nQryEnd) if res.contents.nScore > 30 else None


def main(bam_filename, output_prefix, barcode_seq, cell_barcode_tag, umi_length, new_umi_tag, existing_umi_tag,
         num_read_bases_to_align, ssw_path):
    # Set up our SSW objects:
    ssw = ssw_lib.CSsw(ssw_path)
    alphabet, letter_to_int, mat = ssw_build_matrix()

    # silence message about the .bai file not being found
    pysam.set_verbosity(0)

    # Track the barcode counts:
    out_file_name = f"{output_prefix}.bam"

    num_reads_without_new_umi = 0
    num_reads = 0
    num_new_umis = 0
    num_new_umis_same_as_old = 0

    with pysam.AlignmentFile(bam_filename, 'rb', check_sq=False, require_index=False) as bam_file, \
            tqdm(desc=f"Processing reads", unit="read") as pbar:
        with pysam.AlignmentFile(out_file_name, 'wb', header=bam_file.header) as out_bam_file:

            for read in bam_file:
                pbar.update(1)

                num_reads += 1

                # Let's just align the CBC and UMI
                alignments = get_alignment(ssw,
                                           read.query_sequence[:num_read_bases_to_align],
                                           barcode_seq + read.get_tag(cell_barcode_tag),
                                           alphabet,
                                           letter_to_int,
                                           mat)

                if not alignments:
                    print(f"Could not find adapter + CBC for read: {read.query_name}", file=sys.stderr)
                    num_reads_without_new_umi += 1
                    continue

                # OK, we have our alignment:
                alignment_end = alignments[1]

                new_umi_seq = read.query_sequence[alignment_end+1:alignment_end + 1 + umi_length]
                num_new_umis += 1

                if read.has_tag(existing_umi_tag) and read.get_tag(existing_umi_tag) == new_umi_seq:
                    num_new_umis_same_as_old += 1

                read.set_tag(new_umi_tag, new_umi_seq)

                out_bam_file.write(read)
                pbar.update(1)

    print("Stats:")
    print(f"Total Num Reads: {num_reads}")
    print(f"Num reads with new UMIs: {num_new_umis}")
    print(f"Num reads without new UMIs: {num_reads_without_new_umi}")
    print(f"Num new UMIs same as old UMIs: {num_new_umis_same_as_old}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Align 10x adapter + cell barcode to front of read, get next 10 bases as putative barcode.",
        epilog="")
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument('-b', '--bam', help='BAM filename', required=True)
    requiredNamed.add_argument('-o', "--output-prefix", help='Output prefix', required=True)

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
    parser.add_argument(
        '--num-read-bases-to-align',
        help='Number of bases from the start of each read in which to find the adapter + barcode sequence. .  (default: %(default)d)', type=int,
        default=NUM_BASES_TO_ALIGN
    )

    args = parser.parse_args()

    main(args.bam, args.output_prefix, args.barcode_seq, args.cell_barcode_tag, args.umi_length, args.new_umi_tag,
         args.existing_umi_tag, args.num_read_bases_to_align, args.ssw_path)
