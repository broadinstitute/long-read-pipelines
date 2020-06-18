import re
import argparse
import sys
import pysam

# bamfile1 and bamfile2 MUST originate from the same reads.fastq file
"""
pipeline input
- original_reads.fastq
- ref_asm1
- ref_asm2

1. minimap2(asm1, original_read) -> asm1.bam
2. minimap2(asm2, original_read) -> asm2.bam

asm1.bam and asm2.bam are inputs to this program

1. soft_clipper(original_reads, ref_asm1) -> asm1_softclipped.fastq
2. soft_clipper(orginal_reads ,ref_asm2) -> asm2_softclipped.fastq

 bamfile1 - used ref1 to split reads
 bamfile2 - used ref2 to split reads
"""

def is_primary(flag):
    return (flag & (1 << 8)) + (flag & (1 << 11)) == 0


soft_clips = {}

bam1_unmapped_reads = []
bam = pysam.AlignmentFile(sys.argv[1], "rb")
for alignment in bam:
    if not is_primary(alignment.flag):
        continue
    if alignment.is_unmapped:
        bam1_unmapped_reads.append(alignment.query_name)
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

bam2_unmapped_reads = []
bam = pysam.AlignmentFile(sys.argv[2], "rb")
for alignment in bam:
    if not is_primary(alignment.flag):
        continue
    if alignment.is_unmapped:
        bam2_unmapped_reads.append(alignment.query_name)
        continue

    cigar = alignment.cigarstring

    left_clip_match = re.search('^(\d+)S', cigar)
    left_clip_len = 0 if left_clip_match is None else int(left_clip_match.groups()[0])
    right_clip_match = re.search('(\d+)S$', cigar)
    right_clip_len = 0 if right_clip_match is None else int(right_clip_match.groups()[0])

    # Make sure reads are aligned in the same direction
    if alignment.get_forward_sequence() != alignment.seq:
        left_clip_len, right_clip_len = right_clip_len, left_clip_len

    if alignment.query_name in soft_clips:
        #read also aligned in bam1
        #print(f"{alignment.query_name}: {soft_clips[alignment.query_name]} vs ({left_clip_len}, {right_clip_len})")
        if left_clip_len > 100 and abs(soft_clips[alignment.query_name][0] - left_clip_len) > 100:
            if soft_clips[alignment.query_name][0] > left_clip_len:
                #reference wants to cut a larger chunk
                print(f"{alignment.query_name} (Left) (CutMore) {left_clip_len} vs {soft_clips[alignment.query_name][0]}")
            else:
                print(f"{alignment.query_name} (Left) (CutLess) {left_clip_len} vs {soft_clips[alignment.query_name][0]}")
            # If we were going to clip this side AND the reference would have cut it differently by more than a 100
            # print(f"{alignment.query_name} (Left) {left_clip_len} vs {soft_clips[alignment.query_name][0]}")

        if right_clip_len > 100 and abs(soft_clips[alignment.query_name][1] - right_clip_len) > 100:
            if soft_clips[alignment.query_name][1] > right_clip_len:
                print(f"{alignment.query_name} (Right) (CutMore) {right_clip_len} vs {soft_clips[alignment.query_name][1]}")
            else:
                print(f"{alignment.query_name} (Right) (CutLess) {right_clip_len} vs {soft_clips[alignment.query_name][1]}")

            # If we were going to clip this side AND the reference would have cut it differently by more than a 100
            # print(f"{alignment.query_name} (Right) {right_clip_len} vs {soft_clips[alignment.query_name][1]}")
    elif alignment.query_name in bam1_unmapped_reads:
        if left_clip_len > 100:
            # If we were going to clip this side AND the reference would have cut it differently by more than a 100
            print(f"{alignment.query_name} (Left) (CutLess) {left_clip_len} vs unmapped")
        if right_clip_len > 100:
            # If we were going to clip this side AND the reference would have cut it differently by more than a 100
            print(f"{alignment.query_name} (Right) (CutLess) {right_clip_len} vs unmapped")
        #print(f"{alignment.query_name}: unmapped vs ({left_clip_len}, {right_clip_len})")
    #else:
        #print("ERROR: {alignment.query_name} in {sys.argv[2]} is not in {sys.argv[1]}")

#for unmapped_read in bam2_unmapped_reads:
#    if unmapped_read in soft_clips:
#        print(f"{unmapped_read}: {soft_clips[unmapped_read]} vs unmapped")
#    elif unmapped_read in bam1_unmapped_reads: 
#        print(f"{unmapped_read}: unmapped vs unmapped")
#    else:
#        print("ERROR: {alignment.query_name} in {sys.argv[2]} is not in {sys.argv[1]}")
