import re
import sys
import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--clipping-threshold', dest='clipping_threshold', type=int,
                    help='Minimum soft clip length for read splitting')
parser.add_argument('--ref', dest='ref', type=str,
                    help='Assisting reference')
parser.add_argument('--ref-diff-threshold', dest='ref_diff_threshold', type=int,
                    help='Minimum clipping distance from reference alignment to ignore clipping')
parser.add_argument('--write-ref-conflicts-prefix', dest='write_ref_conflicts_prefix', type=str,
                    help='If set, the tool will only output alignments that are in conflict with the assisting reference '
                         'All lines will be printed in SAM format. The first line will be the unassisted alignment while '
                         'the second line is the assisted alignment from the reference. If you would like to separate the '
                         'unasssisted from the assisted alignments, save the output to a file and pipe it through sed '
                         'to split the file. AN EXAMPLE COMMAND WILL BE ADDED HERE')
args = parser.parse_args()

CLIPPING_THRESHOLD = args.clipping_threshold
REF_DIFF_THRESHOLD = args.ref_diff_threshold
using_ref = args.ref is not None
write_ref_conflicts_prefix = args.write_ref_conflicts_prefix
write_ref_conflicts = write_ref_conflicts_prefix is not None

def is_primary(flag):
    return (flag & (1 << 8)) + (flag & (1 << 11)) == 0

if using_ref:
    soft_clips = {}
    ref_unmapped_reads = []
    ref_alignments_file = pysam.AlignmentFile(args.ref, "r")
    ref_alignments = {}

    for alignment in ref_alignments_file:
        if not is_primary(alignment.flag):
            continue
        if alignment.is_unmapped:
            ref_unmapped_reads.append(alignment.query_name)
            continue


        cigar = alignment.cigarstring
        left_clip_match = re.search('^(\d+)S', cigar)
        left_clip_len = 0 if left_clip_match is None else int(left_clip_match.groups()[0])
        right_clip_match = re.search('(\d+)S$', cigar)
        right_clip_len = 0 if right_clip_match is None else int(right_clip_match.groups()[0])

        # Make sure reads are aligned in the same direction
        if alignment.get_forward_sequence() != alignment.seq:
            left_clip_len, right_clip_len = right_clip_len, left_clip_len

        soft_clips[alignment.query_name] = (left_clip_len, right_clip_len)
        ref_alignments[alignment.query_name] = alignment

    if write_ref_conflicts:
        unassisted_bam = pysam.AlignmentFile(f"{write_ref_conflicts_prefix}_unassisted.bam", "wb", template=ref_alignments_file)
        ref_bam = pysam.AlignmentFile(f"{write_ref_conflicts_prefix}_ref.bam", "wb", template=ref_alignments_file)
    

for line in pysam.AlignmentFile("-", "r"):
    read_id = line.query_name
    cigar = line.cigarstring
    seq = line.seq
    qual = line.qual
    is_reverse = line.is_reverse

    left_clip_len = 0
    right_clip_len = 0
    if cigar is not None:
        left_clip_match = re.search('^(\d+)S', cigar)
        left_clip_len = 0 if left_clip_match is None else int(left_clip_match.groups()[0])
        right_clip_match = re.search('(\d+)S$', cigar)
        right_clip_len = 0 if right_clip_match is None else int(right_clip_match.groups()[0])

    seq_len = len(seq)
    new_primary_start = 0
    new_primary_end = seq_len

    if using_ref:
        if read_id in ref_unmapped_reads:
            ref_left_clip_len = 0
            ref_right_clip_len = 0
        else:
            if is_reverse:
                ref_left_clip_len = soft_clips[read_id][1]
                ref_right_clip_len = soft_clips[read_id][0]
            else:
                ref_left_clip_len = soft_clips[read_id][0]
                ref_right_clip_len = soft_clips[read_id][1]

    if left_clip_len >= CLIPPING_THRESHOLD:
        if using_ref and abs(ref_left_clip_len - left_clip_len) >= REF_DIFF_THRESHOLD:
            #print(f"Skipping due to conflict with aid alignment. Reference wants to cut less")
            if write_ref_conflicts:
                unassisted_bam.write(line)
                ref_bam.write(ref_alignments[read_id])
        else:
            if not write_ref_conflicts:
                print(f"@{read_id}_l\n"
                      f"{seq[:left_clip_len]}\n"
                      f"+\n"
                      f"{qual[:left_clip_len]}")
            new_primary_start = left_clip_len

    if right_clip_len >= CLIPPING_THRESHOLD:
        if using_ref and abs(ref_right_clip_len - right_clip_len) >= REF_DIFF_THRESHOLD:
            #print(f"Skipping due to conflict with aid alignment. Reference wants to cut less")
            if write_ref_conflicts:
                unassisted_bam.write(line)
                ref_bam.write(ref_alignments[read_id])
        else:
            if not write_ref_conflicts:
                print(f"@{read_id}_r\n"
                      f"{seq[seq_len - right_clip_len:]}\n"
                      f"+\n"
                      f"{qual[seq_len - right_clip_len:]}")
            new_primary_end = seq_len - right_clip_len

    if not write_ref_conflicts:
        print(f"@{read_id}\n"
              f"{seq[new_primary_start:new_primary_end]}\n"
              f"+\n"
              f"{qual[new_primary_start:new_primary_end]}")

