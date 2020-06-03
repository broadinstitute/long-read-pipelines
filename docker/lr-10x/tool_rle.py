
import pysam
from Bio import pairwise2
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn3
from Bio.Seq import Seq
from multiprocessing import Pool
import http.client
import time

#bam_filename = '/Users/mgatzen/files/bam/10x_tool/SM-HLD5O/all.corrected.bam'
bam_filename = '/Users/mgatzen/files/bam/10x_tool/ont/1.consensus.bam'

def align(sequence):
    return pairwise2.align.globalms(sequence, ' CTACACGACGCTCTTCCGATCT ', 2, -4, -4, -2, penalize_end_gaps=False)

def position_in_sequence(sequence, pat):
    d = 4
    q = 101
    M = len(pat)
    N = len(sequence)
    i = 0
    j = 0
    p = 0  # hash value for pattern
    t = 0  # hash value for txt
    h = 1

    positions = []

    # The value of h would be "pow(d, M-1)% q"
    for i in range(M - 1):
        h = (h * d) % q

        # Calculate the hash value of pattern and first window
    # of text
    for i in range(M):
        p = (d * p + ord(pat[i])) % q
        t = (d * t + ord(sequence[i])) % q

        # Slide the pattern over text one by one
    for i in range(N - M + 1):
        # Check the hash values of current window of text and
        # pattern if the hash values match then only check
        # for characters on by one
        if p == t:
            # Check for characters one by one
            for j in range(M):
                if sequence[i + j] != pat[j]:
                    break

            j += 1
            # if p == t and pat[0...M-1] = txt[i, i + 1, ...i + M-1]
            if j == M:
                positions.append(i)

        # Calculate hash value for next window of text: Remove
        # leading digit, add trailing digit
        if i < N - M:
            t = (d * (t - ord(sequence[i]) * h) + ord(sequence[i + M])) % q

            # We might get negative values of t, converting it to
            # positive
            if t < 0:
                t = t + q

    return positions

def process_contig_poly_a(contig):

    bam_file = pysam.AlignmentFile(bam_filename, 'rb')

    reads_seen = 0

    max_reads = 10000

    time_start = time.time()
    time_last_segment = time_start

    forward_alignments_histogram, reverse_alignments_histogram = dict(), dict()

    for read in bam_file.fetch(contig=contig):

        if reads_seen >= max_reads:
            break

        reads_seen += 1


        if reads_seen % 1000 == 0:
            print('Processed reads: {}. Last position: {}:{}. Elapsed time: {}'.format(reads_seen, read.reference_name,
                                                                                       read.pos,
                                                                                       time.time() - time_last_segment))
            time_last_segment = time.time()

        forward_sequence = read.seq
        reverse_sequence = str(Seq(forward_sequence).reverse_complement())

        forward_alignments = len(position_in_sequence(forward_sequence, 'TTTTTTTTTTTTTTT'))
        reverse_alignments = len(position_in_sequence(reverse_sequence, 'TTTTTTTTTTTTTTT'))

        if forward_alignments not in forward_alignments_histogram:
            forward_alignments_histogram[forward_alignments] = 0
        forward_alignments_histogram[forward_alignments] += 1

        if reverse_alignments not in reverse_alignments_histogram:
            reverse_alignments_histogram[reverse_alignments] = 0
        reverse_alignments_histogram[reverse_alignments] += 1

    print(forward_alignments_histogram)
    print(reverse_alignments_histogram)


def process_contig(contig, include_odd_sequences):
    bam_file = pysam.AlignmentFile(bam_filename, 'rb')

    reads_seen = 0

    adapter_found, adapter_not_found = 0, 0
    multi_alignments = 0

    forward_adapter_found, forward_adapter_not_found = 0, 0
    forward_multi_alignments = 0
    forward_position_histogram = dict()
    forward_relative_position_histogram = dict()

    reverse_adapter_found, reverse_adapter_not_found = 0, 0
    reverse_multi_alignments = 0
    reverse_position_histogram = dict()
    reverse_relative_position_histogram = dict()

    odd_position_reads = []

    both_ends_aligned = 0

    max_reads = 100000000000

    time_start = time.time()
    time_last_segment = time_start

    for read in bam_file.fetch(contig=contig):

        if reads_seen >= max_reads:
            break

        reads_seen += 1

        if reads_seen % 1000 == 0:
            print('Processed reads: {}. Last position: {}:{}. Elapsed time: {}'.format(reads_seen, read.reference_name, read.pos, time.time() - time_last_segment))
            time_last_segment = time.time()

        forward_sequence = read.seq
        reverse_sequence = str(Seq(forward_sequence).reverse_complement())

        forward_alignments_rk = position_in_sequence(forward_sequence, 'CTACACGACGCTCTTCCGATCT')
        reverse_alignments_rk = position_in_sequence(reverse_sequence, 'CTACACGACGCTCTTCCGATCT')

        for alignment in forward_alignments_rk:
            if alignment not in forward_position_histogram:
                forward_position_histogram[alignment] = 0
            forward_position_histogram[alignment] += 1

            relative_pos = alignment / len(read.seq)
            if relative_pos not in forward_relative_position_histogram:
                forward_relative_position_histogram[relative_pos] = 0
            forward_relative_position_histogram[relative_pos] += 1

            # if alignment > 5:
            #     odd_position_reads.append((read.query_name, read.reference_name, read.pos, forward_alignments_rk, len(forward_sequence), read.seq if include_odd_sequences else None))

        for alignment in reverse_alignments_rk:
            if alignment not in reverse_position_histogram:
                reverse_position_histogram[alignment] = 0
            reverse_position_histogram[alignment] += 1

            relative_pos = alignment / len(read.seq)
            if relative_pos not in reverse_relative_position_histogram:
                reverse_relative_position_histogram[relative_pos] = 0
            reverse_relative_position_histogram[relative_pos] += 1

            # if alignment > 5:
            #     odd_position_reads.append((read.query_name, read.reference_name, read.pos, reverse_alignments_rk, len(reverse_sequence), read.seq if include_odd_sequences else None))

        if len(forward_alignments_rk) == 0:
            forward_adapter_not_found += 1
        else:
            forward_adapter_found += 1
        if len(reverse_alignments_rk) == 0:
            reverse_adapter_not_found += 1
        else:
            reverse_adapter_found += 1

        if len(forward_alignments_rk) == 0 and len(reverse_alignments_rk) == 0:
            adapter_not_found += 1
            continue

        if len(forward_alignments_rk) > 1:
            forward_multi_alignments += 1
        if len(reverse_alignments_rk) > 1:
            reverse_multi_alignments += 1

        if len(forward_alignments_rk) > 0 and len(reverse_alignments_rk) > 0:
            both_ends_aligned += 1
            continue

        if len(forward_alignments_rk) > 1 or len(reverse_alignments_rk) > 1:
            multi_alignments += 1
            continue

        if len(forward_alignments_rk) > 0:
            forward_aligned = True
        else:
            forward_aligned = False

        adapter_found += 1

    print('Total processing time: {}'.format(time.time() - time_start))

    return reads_seen, adapter_found, adapter_not_found, multi_alignments, both_ends_aligned,\
           forward_adapter_found, forward_adapter_not_found, forward_multi_alignments, forward_position_histogram, forward_relative_position_histogram, \
           reverse_adapter_found, reverse_adapter_not_found, reverse_multi_alignments, reverse_position_histogram, reverse_relative_position_histogram, odd_position_reads

def igv_goto(contig, position, length):
    connection = http.client.HTTPConnection('127.0.0.1', 60151)
    connection.send('goto {}:{}'.format(contig, position).encode('utf-8'))
    request = 'region {} {} {}'.format(contig, position, int(position) + int(length))
    #print(request)
    #connection.send('region {} {} {}'.format(contig, position, int(position) + int(length)).encode('utf-8'))
    connection.close()

def click_through_igv(odd_position_reads):
    class print_color:
        PURPLE = '\033[95m'
        CYAN = '\033[96m'
        DARKCYAN = '\033[36m'
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
        END = '\033[0m'

    print('{} odd positioned reads in total.'.format(len(odd_position_reads)))
    for read in odd_position_reads:
        print('{}:{} ({}) {}'.format(read[1], read[2], read[3], read[0]))
        seq = read[5]
        seq_string = ''
        seq_position = 0
        for alignment in read[3]:
            seq_string += seq[seq_position:alignment] + print_color.RED + seq[alignment:alignment+22] + print_color.END
            seq_position = alignment + 22
        seq_string += seq[seq_position:]
        print('{}'.format(seq_string))
        igv_goto(read[1], read[2], read[4])
        input('Press Enter for next position...')

def main():
    # contigs = ['chr{}'.format(num) for num in range(1, 22)]
    # contigs = ['chr1']
    # pool = Pool(len(contigs))
    # results = pool.map(process_contig, contigs)

    use_igv = False

    results = process_contig('chr1', use_igv)
    #results = process_contig_poly_a('chr1')

    #return

    print(results)

    forward_histogram = results[8]
    forward_data = sorted(forward_histogram.items())
    forward_x, forward_y = zip(*forward_data)
    #forward_data_x = np.arange(0, max(forward_x))
    #forward_data_y = [forward_histogram[x] if x in forward_histogram else np.nan for x in range(max(forward_x))]

    forward_relative_histogram = results[9]
    forward_relative_data = sorted(forward_relative_histogram.items())
    forward_relative_x, forward_relative_y = zip(*forward_relative_data)
    forward_relative_y = [(x if x > 0 else np.nan) for x in forward_relative_y]

    reverse_histogram = results[13]
    reverse_data = sorted(reverse_histogram.items())
    reverse_x, reverse_y = zip(*reverse_data)
    reverse_y = [(x if x > 0 else np.nan) for x in reverse_y]

    reverse_relative_histogram = results[14]
    reverse_relative_data = sorted(reverse_relative_histogram.items())
    reverse_relative_x, reverse_relative_y = zip(*reverse_relative_data)
    reverse_relative_x = [(x if x > 0 else np.nan) for x in reverse_relative_x]

    fig, ax = plt.subplots(2, 2, figsize=(10, 6))

    ax[0, 0].plot(forward_x, forward_y)
    ax[0, 0].scatter(forward_x, forward_y, c='red')
    ax[0, 0].set_title('Adapter position in read (forward)')
    ax[0, 0].set_ylim(-2, 20)
    ax[1, 0].plot(reverse_x, reverse_y)
    ax[1, 0].scatter(reverse_x, reverse_y, c='red')
    ax[1, 0].set_title('Adapter position in read (reverse)')
    ax[1, 0].set_ylim(-2, 20)

    ax[0, 1].plot(forward_relative_x, forward_relative_y)
    ax[0, 1].scatter(forward_relative_x, forward_relative_y, c='red')
    ax[0, 1].set_title('Relative adapter position in read (forward)')
    ax[0, 1].set_ylim(-2, 20)
    ax[1, 1].plot(reverse_relative_x, reverse_relative_y)
    ax[1, 1].scatter(reverse_relative_x, reverse_relative_y, c='red')
    ax[1, 1].set_title('Relative adapter position in read (reverse)')
    ax[1, 1].set_ylim(-2, 20)
    plt.show()

    if use_igv:
        click_through_igv(results[15])

main()