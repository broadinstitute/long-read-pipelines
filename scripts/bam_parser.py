import re
import sys

#The read length of a CIGAR is equal to the sum of all numbers excluding soft clips
def get_cigar_length(s):
    return sum(map(int, filter(None, re.sub(r'\D', " ", re.sub(r'\d+S', "", s)).split(" "))))

def is_positive_direction_alignment(sam_flag):
    return not (sam_flag & (1 << (5 - 1)))

for line in sys.stdin:
    tokens = line.split("\t")
  
    read_id = tokens[0]
    sam_flag = int(tokens[1])
    ref_seq = tokens[2]
    start_locus = tokens[3]
    primary_cigar = tokens[5]
    seq = tokens[9]

    supplementary_alignments = []
    for i in range(11, len(tokens)):
        # Assumption - is there only 1 SA flag per primary alignment?
        if tokens[i].startswith("SA:Z:"):
            for st in filter(None, tokens[i][5:].split(";")):
                s_tokens = st.split(",")
                supp_ref = s_tokens[0]
                direction = s_tokens[2]
                locus = s_tokens[1]
                cigar = s_tokens[3]
                mapq = s_tokens[4]

                # Assumption - is the end of the locus always the length of the sequence + the start of the locus?
                supplementary_alignments.append(f'{supp_ref},{direction},{locus}-{int(locus) + int(get_cigar_length(cigar))},{mapq}')
        elif tokens[i].startswith("NM:i"):
            edit_distance = tokens[i].split(":")[2]

    end_locus = int(start_locus) + get_cigar_length(primary_cigar)

    alignment_direction = '+' if is_positive_direction_alignment(sam_flag) else '-'
    print(f'{read_id}\t{alignment_direction}\t{primary_cigar}\t{edit_distance}\t{ref_seq}\t{start_locus}-{end_locus}\t{";".join(supplementary_alignments)}')
