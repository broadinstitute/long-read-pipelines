import re
import sys

CLIPPING_THRESHOLD = 100

for line in sys.stdin:
    tokens = line.split("\t")
    read_id = tokens[0]
    sam_flag = int(tokens[1])
    ref_seq = tokens[2]
    start_locus = tokens[3]
    mapq = tokens[4]
    cigar = tokens[5]
    seq = tokens[9]
    quality = tokens[10]

    left_clip_match = re.search('^(\d+)S', cigar)
    left_clip_len = 0 if left_clip_match is None else int(left_clip_match.groups()[0])
    right_clip_match = re.search('(\d+)S$', cigar)
    right_clip_len = 0 if right_clip_match is None else int(right_clip_match.groups()[0])

    seq_len = len(seq)
    new_primary_start = 0
    new_primary_end = seq_len

    if left_clip_len >= CLIPPING_THRESHOLD:
        print(f"@{read_id}_l\n"
              f"{seq[:left_clip_len]}\n"
              f"+\n"
              f"{quality[:left_clip_len]}")
        new_primary_start = left_clip_len

    if right_clip_len >= CLIPPING_THRESHOLD:
        print(f"@{read_id}_r\n"
              f"{seq[seq_len - right_clip_len:]}\n"
              f"+\n"
              f"{quality[seq_len - right_clip_len:]}")
        new_primary_end = seq_len - right_clip_len

    print(f"@{read_id}\n"
          f"{seq[new_primary_start:new_primary_end]}\n"
          f"+\n"
          f"{quality[new_primary_start:new_primary_end]}")
