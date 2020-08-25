import sys, math, re
import matplotlib.pyplot as plt

print_plots = len(sys.argv) > 2 and sys.argv[2] == "--plot"
f = open(sys.argv[1])

num_primary_alignments = 0
num_has_chimeras = 0
total_chimeras = 0
num_no_chimeras = 0

# Overlap percentage is the % of the chimeric alignment that lies within the primary alignment
chimeric_overlap_percentages = []
chimeric_gap_distance = []
num_chimeras_aligned_to_different_ref_seq = 0

# Clipping counts
MIN_CLIPPING_LENGTHS = [1, 100, 500, 1000]
no_clipping = {}
one_end_clipping = {}
two_end_clipping = {}

# Supplementary Alignment Counts
ALIGN_BACK_OVERLAP_THRESHOLD_PERCENT = 80
CLOSE_BY_GAP_DISTANCE_THRESHOLD_BP = 10000
# true = align in same direction
# false = aligns in opposite direction
aligns_back_to_primary = []
aligns_close_by = []
aligns_far_away = [] # can split by ones that go to different chromosome

mapq_sum = 0
mapq_denom = 0
edit_distance_sum = 0

# regex expressions to parse and count CIGAR strings
misaligned_bases_regex = r"(\d+)([DIHX])"
misaligned_bp_sum = 0
total_soft_clips_length = 0
total_reads_length = 0
read_lengths = []
primary_alignment_lengths = []
per_read_nucleotide_identity = []

unaligned_count = 0
unaligned_bp_sum = 0

for min_clipping_len in MIN_CLIPPING_LENGTHS:
    no_clipping[min_clipping_len] = 0
    one_end_clipping[min_clipping_len] = 0
    two_end_clipping[min_clipping_len] = 0


def get_chimeric_gaps_by_direction(dirs):
    gap_distances = []
    for gap in chimeric_gap_distance:
        if gap[0] in dirs:
            gap_distances.append(gap[1])

    return gap_distances


def compare_sequence_ranges(a_ref, a_range, b_ref, b_range):
    if a_ref != b_ref:
        return False

    a_min = int(a_range.split("-")[0])
    a_max = int(a_range.split("-")[1])

    b_min = int(b_range.split("-")[0])
    b_max = int(b_range.split("-")[1])

    return min(a_max, b_max) - max(a_min, b_min)


def get_seq_range_length(seq_range):
    return int(seq_range.split("-")[1]) - int(seq_range.split("-")[0])


for line in f.readlines():
    tokens = line.split("\t")
    read_id = tokens[0]
    primary_direction = tokens[1]
    cigar = tokens[2]
    edit_distance = tokens[3]
    ref_seq = tokens[4]
    primary_range = tokens[5]
    primary_mapq = tokens[6]
    primary_seq_len = tokens[7]
    chimeras = list(filter(lambda c: get_seq_range_length(c.split(",")[2]) >= 500, filter(None, tokens[8].strip().split(";"))))

    if cigar == "*" and ref_seq == "*":
        unaligned_count = unaligned_count + 1
        unaligned_bp_sum = unaligned_bp_sum + int(primary_seq_len)
        continue

    num_primary_alignments = num_primary_alignments + 1
    total_reads_length = total_reads_length + int(primary_seq_len)
    read_lengths.append(int(primary_seq_len))

    left_clip_match = re.search('^(\d+)S', cigar)
    left_clip_len = 0 if left_clip_match is None else int(left_clip_match.groups()[0])
    right_clip_match = re.search('(\d+)S$', cigar)
    right_clip_len = 0 if right_clip_match is None else int(right_clip_match.groups()[0])
    total_soft_clips_length = total_soft_clips_length + left_clip_len + right_clip_len
    primary_alignment_lengths.append(int(primary_seq_len) - left_clip_len - right_clip_len)

    for MIN_CLIPPING_LENGTH in MIN_CLIPPING_LENGTHS:
        if left_clip_len < MIN_CLIPPING_LENGTH and right_clip_len < MIN_CLIPPING_LENGTH:
            no_clipping[MIN_CLIPPING_LENGTH] = no_clipping[MIN_CLIPPING_LENGTH] + 1
        elif left_clip_len >= MIN_CLIPPING_LENGTH and right_clip_len >= MIN_CLIPPING_LENGTH:
            two_end_clipping[MIN_CLIPPING_LENGTH] = two_end_clipping[MIN_CLIPPING_LENGTH] + 1
        else:
            one_end_clipping[MIN_CLIPPING_LENGTH] = one_end_clipping[MIN_CLIPPING_LENGTH] + 1

    if len(chimeras) > 0: 
        num_has_chimeras = num_has_chimeras + 1
    else:
        num_no_chimeras = num_no_chimeras + 1

    if primary_mapq != 255:
        mapq_sum = mapq_sum + int(primary_mapq)
        mapq_denom = mapq_denom + 1

    edit_distance_sum = edit_distance_sum + int(edit_distance)
    misaligned_bp = sum(int(c) for c, b in re.findall(misaligned_bases_regex, cigar))
    misaligned_bp_sum = misaligned_bp_sum + misaligned_bp

    read_len_without_clipping = int(primary_seq_len) - left_clip_len - right_clip_len
    nucleotide_identity = (read_len_without_clipping - misaligned_bp) * 100 / read_len_without_clipping
    per_read_nucleotide_identity.append(nucleotide_identity)

    for chimera in chimeras:
        total_chimeras = total_chimeras + 1

        chimera_tokens = chimera.split(",")
        chimera_ref = chimera_tokens[0]
        chimera_dir = chimera_tokens[1]
        chimera_range = chimera_tokens[2]

        overlap = compare_sequence_ranges(chimera_ref, chimera_range, ref_seq, primary_range)

        if overlap is not False:
            if overlap > 0:
                chimeric_overlap_percentages.append(overlap * 100 / get_seq_range_length(chimera_range))
            else:
                chimeric_gap_distance.append((chimera_dir, overlap * -1))
        else:
            num_chimeras_aligned_to_different_ref_seq = num_chimeras_aligned_to_different_ref_seq + 1

        if overlap is not False:
            if overlap > 0:
                overlap_percent = overlap * 100 / get_seq_range_length(chimera_range)

                if overlap_percent >= ALIGN_BACK_OVERLAP_THRESHOLD_PERCENT:
                    aligns_back_to_primary.append(chimera_dir == primary_direction)
                else:
                    aligns_close_by.append(chimera_dir == primary_direction)
            else:
                if overlap * -1 <= CLOSE_BY_GAP_DISTANCE_THRESHOLD_BP:
                    aligns_close_by.append(chimera_dir == primary_direction)
                else:
                    aligns_far_away.append(chimera_dir == primary_direction)
        else:
            aligns_far_away.append(chimera_dir == primary_direction)


def percent_string(num, denom):
    return f"{format(num * 100 / denom, '.2f')}% ({num})"

print("\n\n\n")
print(f"Number of reads: {num_primary_alignments + unaligned_count}")
print(f"Aligned reads: {percent_string(num_primary_alignments, num_primary_alignments + unaligned_count)}")
print(f"Unaligned reads: {percent_string(unaligned_count, unaligned_count + num_primary_alignments)}")
print(f"Average Read length (including unaligned): {math.floor((total_reads_length + unaligned_bp_sum)/(num_primary_alignments + unaligned_count))}")
print(f"Average Read length (aligned only): {math.floor(total_reads_length / num_primary_alignments)}")
print(f"Average Primary alignment length (no soft clips): {math.floor((total_reads_length - total_soft_clips_length) / num_primary_alignments)}")

print("Clipping Summary")
print(f"Clipping Threshold\tNo Clipping\t\tOne Sided Clip\t\tTwo Sided Clips")
for MIN_CLIPPING_LENGTH in MIN_CLIPPING_LENGTHS:
    print(f"{MIN_CLIPPING_LENGTH}"
          f"\t\t\t{percent_string(no_clipping[MIN_CLIPPING_LENGTH], num_primary_alignments)}"
          f"\t\t{percent_string(one_end_clipping[MIN_CLIPPING_LENGTH], num_primary_alignments)}"
          f"\t\t{percent_string(two_end_clipping[MIN_CLIPPING_LENGTH], num_primary_alignments)}")

print("\n")

print("Supplementary Reads Summary")
print(f"Reads with supplementary alignments: {percent_string(num_has_chimeras, num_primary_alignments)}")
print(f"Aligns back to primary: {percent_string(len(aligns_back_to_primary), total_chimeras)}")
print(f"   in same direction: "
      f"{percent_string(aligns_back_to_primary.count(True), len(aligns_back_to_primary))}")
print(f"   in reverse direction: "
      f"{percent_string(aligns_back_to_primary.count(False), len(aligns_back_to_primary))}")
print("")


print(f"Aligns close by: {percent_string(len(aligns_close_by), total_chimeras)}")
print(f"   in same direction: "
      f"{percent_string(aligns_close_by.count(True), len(aligns_close_by))}")
print(f"   in reverse direction: "
      f"{percent_string(aligns_close_by.count(False), len(aligns_close_by))}")
print("")


print(f"Aligns far away: {percent_string(len(aligns_far_away), total_chimeras)}")
print(f"   in same direction: "
      f"{percent_string(aligns_far_away.count(True), len(aligns_far_away))}")
print(f"   in reverse direction: "
      f"{percent_string(aligns_far_away.count(False), len(aligns_far_away))}")
print("")

print("\n")
print(f"Average MAPQ: {mapq_sum / mapq_denom}")
total_bp_without_soft_clips = total_reads_length - total_soft_clips_length
print(f"Average Nucleotide Identity* of primary alignments: {percent_string(total_bp_without_soft_clips - misaligned_bp_sum, total_bp_without_soft_clips)}")
print(f"Average Nucleotide Identity* of primary alignments with soft clips: {percent_string(total_reads_length - misaligned_bp_sum - total_soft_clips_length, total_reads_length)}")
print(f">80%: {percent_string(sum(p >= 80 for p in per_read_nucleotide_identity), num_primary_alignments)}")
print(f">85%: {percent_string(sum(p >= 85 for p in per_read_nucleotide_identity), num_primary_alignments)}")
print(f">90%: {percent_string(sum(p >= 90 for p in per_read_nucleotide_identity), num_primary_alignments)}")
print(f">95%: {percent_string(sum(p >= 95 for p in per_read_nucleotide_identity), num_primary_alignments)}")
print("* edit distance to ref / read length")

plt.hist(read_lengths, 50)
plt.title("Read Lengths")
plt.show()

#fig, axs = plt.subplots(2)
#axs[0].hist(read_lengths, 50)
#axs[1].hist(primary_alignment_lengths, 50)
#plt.show()
