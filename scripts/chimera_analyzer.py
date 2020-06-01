import sys
import statistics
import plotille
import matplotlib.pyplot as plt
import re

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
CLOSE_BY_GAP_DISTANCE_THRESHOLD_BP = 1000
# true = align in same direction
# false = aligns in opposite direction
aligns_back_to_primary = []
aligns_close_by = []
aligns_far_away = [] # can split by ones that go to different chromosome

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
    chimeras = list(filter(lambda c: get_seq_range_length(c.split(",")[2]) >= 500, filter(None, tokens[6].strip().split(";"))))

    left_clip_match = re.search('^(\d+)S', cigar)
    left_clip_len = 0 if left_clip_match is None else int(left_clip_match.groups()[0])
    right_clip_match = re.search('(\d+)S$', cigar)
    right_clip_len = 0 if right_clip_match is None else int(right_clip_match.groups()[0])

    for MIN_CLIPPING_LENGTH in MIN_CLIPPING_LENGTHS:
        if left_clip_len < MIN_CLIPPING_LENGTH and right_clip_len < MIN_CLIPPING_LENGTH:
            no_clipping[MIN_CLIPPING_LENGTH] = no_clipping[MIN_CLIPPING_LENGTH] + 1
        elif left_clip_len >= MIN_CLIPPING_LENGTH and right_clip_len >= MIN_CLIPPING_LENGTH:
            two_end_clipping[MIN_CLIPPING_LENGTH] = two_end_clipping[MIN_CLIPPING_LENGTH] + 1
        else:
            one_end_clipping[MIN_CLIPPING_LENGTH] = one_end_clipping[MIN_CLIPPING_LENGTH] + 1

    num_primary_alignments = num_primary_alignments + 1
    if len(chimeras) > 0: 
        num_has_chimeras = num_has_chimeras + 1
    else:
        num_no_chimeras = num_no_chimeras + 1

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

print(f"Primary Alignments: {num_primary_alignments}")
print(f"Alignments with Chimeras: {num_has_chimeras} ({num_has_chimeras * 100 / num_primary_alignments} %)")
print(f"Average # of chimeras out of alignments with at least one: {total_chimeras / num_has_chimeras}")

print("")
print(f"Total number of Chimeras: {total_chimeras}")
print(f"Total Chimeras overlapping onto primary alignment: {len(chimeric_overlap_percentages)} ({len(chimeric_overlap_percentages) * 100 / total_chimeras} %)")
print(f"\tAverage % of overlap. (How much of the chimeric alignment is in the primary): {statistics.mean(chimeric_overlap_percentages)} %")
if print_plots:
    print(plotille.histogram(chimeric_overlap_percentages))

print(f"Total Chimeras with no overlap: {len(chimeric_gap_distance) + num_chimeras_aligned_to_different_ref_seq} ({(len(chimeric_gap_distance) + num_chimeras_aligned_to_different_ref_seq) * 100 / total_chimeras} %)")
print(f"\tTotal Chimeras aligning to same chromosome: {len(chimeric_gap_distance)}")
print(f"\tAverage gap distance between chimeras aligned to same chromosome: {statistics.mean(get_chimeric_gaps_by_direction(['+', '-']))}")

if print_plots:
    positive_plot = plt.figure(1)
    plt.title("Positive Alignments")
    plt.xlim(0, 5000)
    plt.xlabel("Gap Distance (bp)")
    plt.ylabel("Count")
    plt.hist(get_chimeric_gaps_by_direction(['+']), bins="auto")

    negative_plot = plt.figure(2)
    plt.title("Reversed Alignments")
    plt.xlim(0, 5000)
    plt.xlabel("Gap Distance (bp)")
    plt.ylabel("Count")
    plt.hist(get_chimeric_gaps_by_direction(['-']), bins="auto")

    all_plot = plt.figure(3)
    plt.title("All Alignments")
    plt.xlim(0, 5000)
    plt.xlabel("Gap Distance (bp)")
    plt.ylabel("Count")
    all_gaps = get_chimeric_gaps_by_direction(['-', '+'])
    plt.hist(all_gaps, bins=range(min(all_gaps), max(all_gaps) + 100, 100))

    plt.show()

print(f"\tTotal Chimeras aligning to different chromosome: {num_chimeras_aligned_to_different_ref_seq}")

def percent_string(num, denom):
    return format(num * 100 / denom, '.2f') + '%'

print("\n\n\n")
print("Clipping Summary")
print(f"Clipping Threshold\tNo Clipping\t\tOne Sided Clip\t\tTwo Sided Clips")
for MIN_CLIPPING_LENGTH in MIN_CLIPPING_LENGTHS:
    print(f"{MIN_CLIPPING_LENGTH}"
          f"\t\t\t{percent_string(no_clipping[MIN_CLIPPING_LENGTH], num_primary_alignments)} ({no_clipping[MIN_CLIPPING_LENGTH]})"
          f"\t\t{percent_string(one_end_clipping[MIN_CLIPPING_LENGTH], num_primary_alignments)} ({one_end_clipping[MIN_CLIPPING_LENGTH]}) "
          f"\t\t{percent_string(two_end_clipping[MIN_CLIPPING_LENGTH], num_primary_alignments)} ({two_end_clipping[MIN_CLIPPING_LENGTH]})")

print("\n")

print(f"Aligns back to primary: {percent_string(len(aligns_back_to_primary), total_chimeras)} ({len(aligns_back_to_primary)})")
print(f"   in same direction: "
      f"{percent_string(aligns_back_to_primary.count(True), len(aligns_back_to_primary))} ({aligns_back_to_primary.count(True)})")
print(f"   in reverse direction: "
      f"{percent_string(aligns_back_to_primary.count(False), len(aligns_back_to_primary))} ({aligns_back_to_primary.count(False)})")
print("")

print(f"Aligns close by: {percent_string(len(aligns_close_by), total_chimeras)} ({len(aligns_close_by)})")
print(f"   in same direction: "
      f"{percent_string(aligns_close_by.count(True), len(aligns_close_by))} ({aligns_close_by.count(True)})")
print(f"   in reverse direction: "
      f"{percent_string(aligns_close_by.count(False), len(aligns_close_by))} ({aligns_close_by.count(False)})")
print("")

print(f"Aligns far away: {percent_string(len(aligns_far_away), total_chimeras)} ({len(aligns_far_away)})")
print(f"   in same direction: "
      f"{percent_string(aligns_far_away.count(True), len(aligns_far_away))} ({aligns_far_away.count(True)})")
print(f"   in reverse direction: "
      f"{percent_string(aligns_far_away.count(False), len(aligns_far_away))} ({aligns_far_away.count(False)})")
print("")
