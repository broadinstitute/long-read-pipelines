import re
import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--clipping-threshold', dest='clipping_threshold', type=int,
                    help='Minimum soft clip length for read splitting')
parser.add_argument('--ref', dest='ref', type=str,
                    help='Assisting reference')
args = parser.parse_args()

CLIPPING_THRESHOLD = args.clipping_threshold
using_ref = args.ref is not None

def is_primary(flag):
    return (flag & (1 << 8)) + (flag & (1 << 11)) == 0

# Stores (left_clip, right_clip) for the forward strand

if using_ref:
    soft_clips = {}
    ref_unmapped_reads = []
    ref_alignments = pysam.AlignmentFile(args.ref, "rb")
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

for alignment in pysam.AlignmentFile("-", "r"):
    if not is_primary(alignment.flag):
        continue

    cigar = alignment.cigarstring

    left_clip_len = 0
    right_clip_len = 0
    if cigar is not None:
        left_clip_match = re.search('^(\d+)S', cigar)
        left_clip_len = 0 if left_clip_match is None else int(left_clip_match.groups()[0])
        right_clip_match = re.search('(\d+)S$', cigar)
        right_clip_len = 0 if right_clip_match is None else int(right_clip_match.groups()[0])

    seq_len = len(alignment.seq)
    new_primary_start = 0
    new_primary_end = seq_len

    if alignment.query_name not in soft_clips and alignment.query_name in ref_unmapped_reads:
        ref_left_clip_len = 0
        ref_right_clip_len = 0
    else:
        if alignment.is_reverse:
            ref_left_clip_len = soft_clips[alignment.query_name][1]
            ref_right_clip_len = soft_clips[alignment.query_name][0]
        else:
            ref_left_clip_len = soft_clips[alignment.query_name][0]
            ref_right_clip_len = soft_clips[alignment.query_name][1]

    if left_clip_len >= CLIPPING_THRESHOLD:
        if using_ref and abs(ref_left_clip_len - left_clip_len) >= CLIPPING_THRESHOLD:
            if ref_left_clip_len < left_clip_len:
                #print("Skipping where ref says to clip less")
                pass
            else:
                #print("Skipping where ref says to clip more")
                pass
        else:
            print(f"@{alignment.query_name}_l\n"
                  f"{alignment.seq[:left_clip_len]}\n"
                  f"+\n"
                  f"{alignment.qual[:left_clip_len]}")
            new_primary_start = left_clip_len

    if right_clip_len >= CLIPPING_THRESHOLD:
        if using_ref and abs(ref_right_clip_len - right_clip_len) >= CLIPPING_THRESHOLD:
            if ref_right_clip_len < right_clip_len:
                #print("Skipping where ref says to clip less")
                pass
            else:
                #print("Skipping where ref says to clip more")
                pass
        else:
            print(f"@{alignment.query_name}_r\n"
                  f"{alignment.seq[seq_len - right_clip_len:]}\n"
                  f"+\n"
                  f"{alignment.qual[seq_len - right_clip_len:]}")
            new_primary_end = seq_len - right_clip_len

    print(f"@{alignment.query_name}\n"
          f"{alignment.seq[new_primary_start:new_primary_end]}\n"
          f"+\n"
          f"{alignment.qual[new_primary_start:new_primary_end]}")

