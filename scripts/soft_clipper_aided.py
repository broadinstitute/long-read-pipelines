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
parser.add_argument('--keep-all-ref-conflicts', dest='keep_all_ref_conflicts', action='store_true',
                    help='Still soft clips regions that differ from the aided reference IF the reference wants to cut a greater region')
args = parser.parse_args()

CLIPPING_THRESHOLD = args.clipping_threshold
REF_DIFF_THRESHOLD = args.ref_diff_threshold
keep_all_ref_conflicts = args.keep_all_ref_conflicts
using_ref = args.ref is not None

def is_primary(flag):
    return (flag & (1 << 8)) + (flag & (1 << 11)) == 0

if using_ref:
    soft_clips = {}
    ref_unmapped_reads = []
    ref_alignments = pysam.AlignmentFile(args.ref, "r")

    for alignment in ref_alignments:
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


for line in sys.stdin:
    tokens = line.split("\t")

    read_id = tokens[0]
    sam_flag = int(tokens[1])
    cigar = tokens[5]
    seq = tokens[9]
    qual = tokens[10]
    is_reverse = (sam_flag & (1 << 4)) != 0

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
        if using_ref and abs(ref_left_clip_len - left_clip_len) >= REF_DIFF_THRESHOLD and (ref_left_clip_len < left_clip_len or keep_all_ref_conflicts):
            #print(f"Skipping due to conflict with aid alignment. Reference wants to cut less")
            pass
        else:
            print(f"@{read_id}_l\n"
                  f"{seq[:left_clip_len]}\n"
                  f"+\n"
                  f"{qual[:left_clip_len]}")
            new_primary_start = left_clip_len

    if right_clip_len >= CLIPPING_THRESHOLD:
        if using_ref and abs(ref_right_clip_len - right_clip_len) >= REF_DIFF_THRESHOLD and (ref_right_clip_len < right_clip_len or keep_all_ref_conflicts):
            #print(f"Skipping due to conflict with aid alignment. Reference wants to cut less")
            pass
        else:
            print(f"@{read_id}_r\n"
                  f"{seq[seq_len - right_clip_len:]}\n"
                  f"+\n"
                  f"{qual[seq_len - right_clip_len:]}")
            new_primary_end = seq_len - right_clip_len

    print(f"@{read_id}\n"
          f"{seq[new_primary_start:new_primary_end]}\n"
          f"+\n"
          f"{qual[new_primary_start:new_primary_end]}")

