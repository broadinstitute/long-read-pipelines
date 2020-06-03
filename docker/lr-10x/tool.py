import subprocess

import pysam
import argparse
from Bio.Seq import Seq
import time
import os.path
import gzip
import ctypes
import sys
sys.path.append('/lrma')
from ssw import ssw_lib

BARCODE_LENGTH = 16
UMI_LENGTH = 12
POLY_T_LENGTH = 10

ADAPTER_TAG = 'ZA'
BARCODE_TAG = 'CB'
RAW_BARCODE_TAG = 'CR'
UMI_TAG = 'ZU'

def read_barcodes(barcodes_filename):
    """
    Reads a line-separated file of barcodes. If a line ends with '-1', this suffix is removed. A gzip file is supported if the extension of the file is .gz
    :param barcodes_filename: Filename of the barcode file.
    :return: A set of barcodes.
    """
    barcodes = set()
    _, barcodes_file_extension = os.path.splitext(barcodes_filename)
    is_gzip = barcodes_file_extension == '.gz'
    with (open(barcodes_filename) if not is_gzip else gzip.open(barcodes_filename)) as barcodes_file:
        for line in barcodes_file:
            if is_gzip:
                line = line.decode('utf-8')
            if line[-3:-1] == '-1':
                barcodes.add(line[:-3])
            else:
                barcodes.add(line[:-1])
    return barcodes


def find_barcode(sequence, adapter_alignment_end):
    """
    Infers the barcode and its position given a sequence and the position of the end of the adapter.
    :param sequence: The sequence of the read end
    :param adapter_alignment_end: The position of the end of the adapter alignment
    :return: raw barcode, position of the beginning of the barcode in the sequence
    """
    estimated_barcode_position = adapter_alignment_end + 1
    if estimated_barcode_position + BARCODE_LENGTH >= len(sequence):
        return None, None
    estimated_barcode = sequence[estimated_barcode_position:estimated_barcode_position + BARCODE_LENGTH]
    return estimated_barcode, estimated_barcode_position


def find_umi(sequence, barcode_position):
    """
    Infers the UMI and its position given a sequence and the position of the start of the barcode
    :param sequence: The sequence of the read end
    :param barcode_position: position of the beginning of the barcode sequence
    :return: raw UMI, position of the beginning of the UMI in the sequence
    """
    umi_position = barcode_position + BARCODE_LENGTH
    if umi_position + UMI_LENGTH >= len(sequence):
        return None, None

    umi = sequence[umi_position:umi_position + UMI_LENGTH]

    return umi, umi_position


def check_poly_t_ratio(sequence, umi_position):
    """
    Returns the ratio of Ts in a section of length POLY_T_LENGTH after the UMI in the read end
    :param sequence: The sequence of the read end
    :param umi_position: The position of the beginning of the UMI in the sequence
    :return: The ratio of Ts in a section of length POLY_T_LENGTH after the UMI in the read end
    """
    if umi_position + UMI_LENGTH + POLY_T_LENGTH >= len(sequence):
        return None
    estimated_poly_t = sequence[umi_position + UMI_LENGTH:umi_position + UMI_LENGTH + POLY_T_LENGTH]
    return estimated_poly_t.count('T') / POLY_T_LENGTH



class AnalysisStats:
    """
    Class for storing stats during processing of a file. Each parameter has to be specified in 3 places: the init
    function, in the header string, and in the print_string function.
    """
    def __init__(self):
        self.reads_seen = 0
        self.adapter_found = 0
        self.adapter_not_found = 0
        self.adapter_in_both_ends = 0
        self.barcode_found = 0
        self.barcode_not_corrected = 0
        self.barcode_corrected = 0
        self.barcode_filtered = 0
        self.barcode_not_found = 0
        self.umi_found = 0
        self.umi_not_found = 0
        self.poly_t_found = 0
        self.poly_t_not_found = 0
        self.sequence_too_short_for_umi = 0
        self.sequence_too_short_for_poly_t = 0

        self.forward_tso_in_not_found = 0
        self.reverse_tso_in_not_found = 0
        self.both_tsos_in_not_found = 0
        self.forward_tso_in_forward_found = 0
        self.reverse_tso_in_forward_found = 0
        self.both_tsos_in_forward_found = 0
        self.forward_tso_in_reverse_found = 0
        self.reverse_tso_in_reverse_found = 0
        self.both_tsos_in_reverse_found = 0

        self.barcode_in_10x_and_illumina_whitelist = 0
        self.barcode_in_10x_whitelist = 0
        self.barcode_in_illumina_whitelist = 0
        self.barcode_not_in_any_whitelist = 0

        self.starcode_clusters = 0
        self.corrected_from_no_list_to_10x = 0

    print_header = ('reads_seen\t'
                    'adapter_found\t'
                    'adapter_not_found\t'
                    'adapter_in_both_ends\t'
                    'barcode_found\t'
                    'barcode_not_corrected\t'
                    'barcode_corrected\t'
                    'barcode_filtered\t'
                    'barcode_not_found\t'
                    'umi_found\t'
                    'umi_not_found\t'
                    'poly_t_found\t'
                    'poly_t_not_found\t'
                    'sequence_too_short_for_umi\t'
                    'sequence_too_short_for_poly_t\t'
                    'forward_tso_in_not_found\t'
                    'reverse_tso_in_not_found\t'
                    'both_tsos_in_not_found\t'
                    'forward_tso_in_forward_found\t'
                    'reverse_tso_in_forward_found\t'
                    'both_tsos_in_forward_found\t'
                    'forward_tso_in_reverse_found\t'
                    'reverse_tso_in_reverse_found\t'
                    'both_tsos_in_reverse_found\t'
                    'barcode_in_10x_and_illumina_whitelist\t'
                    'barcode_in_10x_whitelist\t'
                    'barcode_in_illumina_whitelist\t'
                    'barcode_not_in_any_whitelist\t'
                    'starcode_clusters\t'
                    'corrected_from_no_list_to_10x\n')

    def print_string(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            self.reads_seen,
            self.adapter_found,
            self.adapter_not_found,
            self.adapter_in_both_ends,
            self.barcode_found,
            self.barcode_not_corrected,
            self.barcode_corrected,
            self.barcode_filtered,
            self.barcode_not_found,
            self.umi_found,
            self.umi_not_found,
            self.poly_t_found,
            self.poly_t_not_found,
            self.sequence_too_short_for_umi,
            self.sequence_too_short_for_poly_t,
            self.forward_tso_in_not_found,
            self.reverse_tso_in_not_found,
            self.both_tsos_in_not_found,
            self.forward_tso_in_forward_found,
            self.reverse_tso_in_forward_found,
            self.both_tsos_in_forward_found,
            self.forward_tso_in_reverse_found,
            self.reverse_tso_in_reverse_found,
            self.both_tsos_in_reverse_found,
            self.barcode_in_10x_and_illumina_whitelist,
            self.barcode_in_10x_whitelist,
            self.barcode_in_illumina_whitelist,
            self.barcode_not_in_any_whitelist,
            self.starcode_clusters,
            self.corrected_from_no_list_to_10x)


def ssw_build_matrix(match_score=2, mismatch_score=1):
    """
    Builds a matrix of ctypes values for matches and mismatches for the Smith-Waterman algorithm
    :param match_score: The score for matches
    :param mismatch_score: The penalty for mismatches, i.e. a positive number results in a penalty of that number.
    :return: A matrix for use in the ssw algorithm
    """
    alphabet = ['A', 'C', 'G', 'T', 'N']
    dRc = {'A': 'C', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'C', 'c': 'G', 'g': 'C', 't': 'A'}
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
    :return: The position of the last base of the adapter sequence in the read end
    """
    sequence_numbers = to_int(sequence, alphabet, letter_to_int)
    adapter_numbers = to_int(adapter_sequence, alphabet, letter_to_int)

    flag = 1
    mask_length = len(adapter_sequence) // 2 if len(adapter_sequence) >= 30 else 15

    qProfile = ssw.ssw_init(sequence_numbers, ctypes.c_int32(len(sequence)), mat, len(alphabet), 2)

    res = ssw.ssw_align(qProfile, adapter_numbers, ctypes.c_int32(len(adapter_sequence)), open_penalty, extension_penalty, flag, 0, 0, mask_length)

    return res.contents.nQryEnd if res.contents.nScore > 30 else None


def align(read, stats, ssw, alphabet, letter_to_int, mat, read_end_length, adapter_sequence, tso_sequence):
    """
    Performs the alignment of the read end to the adapter sequence and the TSO sequence
    :param read: The read object
    :param stats: The AnalysisStats object
    :param ssw: The ssw object
    :param alphabet: The alphabet object for the ssw algorithm
    :param letter_to_int: The dict to convert base letters to numbers
    :param mat: The match-mismatch score matrix for the ssw algorithm
    :param read_end_length: Length of the read end that must include the adapter, barcode, UMI, and poly-T tail (recommended for PacBio: 80, recommended for Oxford Nanopore: 250)
    :param adapter_sequence: The adapter sequence to align to
    :param tso_sequence: The TSO sequence to align to
    :return: The sequence of the read end that the adapter was found in, the position of the last base of the adapter sequence in the read end
    """
    read_seq = read.seq
    read_seq_reversed = str(Seq(read_seq).reverse_complement())
    five_prime_end = read_seq[:read_end_length]
    three_prime_end_reversed = read_seq_reversed[:read_end_length]

    five_prime_alignment_end = get_alignment(ssw, five_prime_end, adapter_sequence, alphabet, letter_to_int, mat)
    three_prime_alignment_end = get_alignment(ssw, three_prime_end_reversed, adapter_sequence, alphabet, letter_to_int, mat)

    five_prime_tso_alignment_end = get_alignment(ssw, five_prime_end, tso_sequence, alphabet, letter_to_int, mat)
    three_prime_tso_alignment_end = get_alignment(ssw, three_prime_end_reversed, tso_sequence, alphabet, letter_to_int, mat)

    if five_prime_alignment_end is None and three_prime_alignment_end is None:
        stats.adapter_not_found += 1

        if five_prime_tso_alignment_end:
            if three_prime_tso_alignment_end:
                stats.both_tsos_in_not_found += 1
            else:
                stats.forward_tso_in_not_found += 1
        elif three_prime_tso_alignment_end:
            stats.reverse_tso_in_not_found += 1
        return None, None

    if five_prime_alignment_end and three_prime_alignment_end:
        stats.adapter_in_both_ends += 1
        return None, None

    if five_prime_alignment_end:
        alignment_end = five_prime_alignment_end
        sequence = read_seq

        if five_prime_tso_alignment_end:
            if three_prime_tso_alignment_end:
                stats.both_tsos_in_forward_found += 1
            else:
                stats.forward_tso_in_forward_found += 1
        elif three_prime_tso_alignment_end:
            stats.reverse_tso_in_forward_found += 1
    else:
        alignment_end = three_prime_alignment_end
        sequence = read_seq_reversed

        if five_prime_tso_alignment_end:
            if three_prime_tso_alignment_end:
                stats.both_tsos_in_reverse_found += 1
            else:
                stats.forward_tso_in_reverse_found += 1
        elif three_prime_tso_alignment_end:
            stats.reverse_tso_in_reverse_found += 1

    # Valid adapter found
    stats.adapter_found += 1
    return sequence, alignment_end


def process_barcode(sequence, adapter_alignment_end, stats, whitelist_10x, whitelist_illumina):
    """
    Finds the barcode and determines if it is present in the provded whitelists
    :param sequence: The sequence of the read end
    :param adapter_alignment_end: The position of the last base of the adapter sequence in the read end
    :param stats: The AnalysisStats object
    :param whitelist_10x: Whitelist provided by 10x. Can be None, unless whitelist_illumina is None.
    :param whitelist_illumina: Whitelist provided by short read sequencing. Can be None
    :return: The raw barcode sequence, the position of the first read of the barcode in the read end sequence, whether or not the barcode is present in the 10x whitelist, whether or not the barcode is present in the Illumina whitelist
    """
    barcode, barcode_position = find_barcode(sequence, adapter_alignment_end)

    if not barcode:
        stats.barcode_not_found += 1
        return None, None, None, None

    if whitelist_10x is not None and barcode in whitelist_10x:
        barcode_in_10x = True
    else:
        barcode_in_10x = False
    if whitelist_illumina is not None and barcode in whitelist_illumina:
        barcode_in_illumina = True
    else:
        barcode_in_illumina = False

    if barcode_in_10x:
        if barcode_in_illumina:
            stats.barcode_in_10x_and_illumina_whitelist += 1
        else:
            stats.barcode_in_10x_whitelist += 1
    else:
        if barcode_in_illumina:
            stats.barcode_in_illumina_whitelist += 1
        else:
            stats.barcode_not_in_any_whitelist += 1

    stats.barcode_found += 1

    return barcode, barcode_position, barcode_in_10x, barcode_in_illumina


def process_umi(sequence, barcode_position, stats):
    """
    Infers the UMI and its position
    :param sequence: The sequence of the read end
    :param barcode_position: The position of the first base of the barcode in the read end
    :param stats: The AnalysisStats object
    :return: The raw UMI sequence, the position of the first base of the UMI in the read end sequence
    """
    umi, umi_position = find_umi(sequence, barcode_position)
    if not umi:
        stats.sequence_too_short_for_umi += 1
        return None, None

    stats.umi_found += 1

    return umi, umi_position


def process_poly_t(sequence, umi_position, stats):
    """
    Checks if the poly-T tail is present by testing if the ratio of Ts in the section of a certain length in the read end after the UMI is greater than a given number.
    :param sequence: The sequence of the read end
    :param umi_position: The position of the first base of the UMI in the read end
    :param stats: The AnalysisStats object
    :return: True if a poly-T tail is found
    """
    poly_t_ratio = check_poly_t_ratio(sequence, umi_position)
    if not poly_t_ratio:
        stats.sequence_too_short_for_poly_t += 1
        return False
    elif poly_t_ratio < 0.9:
        stats.poly_t_not_found += 1
        return False
    else:
        stats.poly_t_found += 1
        return True


def perform_barcode_correction_starcode(observed_barcodes, analysis_name, starcode_path, stats=None):
    """
    Performs the barcode correction using starcode
    :param observed_barcodes: dict of observed barcodes with the barcode sequences as keys and the number of observations as values
    :param analysis_name: Prefix for the stats files
    :param starcode_path: Relative or absolute path to the starcode executable
    :param stats: The AnalysisStats file
    :return: A dictionary with the raw barcodes as keys and the corresponding corrected barcodes as values
    """
    correction_dict = dict()

    starcode_args = [starcode_path, '--print-clusters', '--quiet']

    starcode_proc = subprocess.Popen(starcode_args, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    encoding = 'ascii'

    for (barcode, occurrence) in observed_barcodes.items():
        for i in range(occurrence):
            starcode_proc.stdin.write(bytearray(barcode + '\n', encoding=encoding))

    # Run starcode
    starcode_proc.stdin.close()

    starcode_output_file = open(analysis_name + '_starcode.tsv', 'w') if analysis_name is not None else None
    if starcode_output_file is not None:
        starcode_output_file.write('cluster\toccurrence\tindices\n')

    num_clusters = 0
    for cluster_line in starcode_proc.stdout:
        num_clusters += 1
        cluster_line = cluster_line.decode(encoding)
        cluster_data = cluster_line.strip().split('\t')
        cluster = cluster_data[0]
        cluster_sequences = cluster_data[2].split(',')
        for cluster_sequence in cluster_sequences:
            correction_dict[cluster_sequence] = cluster
        if starcode_output_file is not None:
            starcode_output_file.write(cluster_line)
    if starcode_output_file is not None:
        starcode_output_file.close()

    if stats is not None:
        stats.starcode_clusters = num_clusters

    return correction_dict


def main(bam_filename, analysis_name, adapter_fasta_filename, tso_fasta_filename, whitelist_10x_filename, whitelist_illumina_filename, max_reads, contig, read_end_length, record_umis, ssw_path, starcode_path):
    """
    Main function for the tool
    :param bam_filename: Filename of the reads BAM file
    :param analysis_name: Prefix for storing the analysis stats files
    :param adapter_fasta_filename: Filename for the FASTA file of the adapter sequence
    :param tso_fasta_filename: Filename for the FASTA file of the TSO sequence
    :param whitelist_10x_filename: Filename of the 10x whitelist. Can be None or empty, unless whitelist_illumina_filename is None or empty
    :param whitelist_illumina_filename: Filename of the Illumina whitelist. Can be None or empty
    :param max_reads: Number of reads to process before stopping the processing. Can be None to process the entire file.
    :param contig: Contig name to limit the processing to. Can be None to process all contigs.
    :param read_end_length: Number of bases to look for the adapter, barcode, UMI, and poly-T tail in both ends of each read. Recommendatios for {PacBio: 80, Oxford Nanopore: 250)
    :param record_umis: Whether or not to store the UMIs in a separate file
    :param ssw_path: Path to the ssw library
    :param starcode_path: Path to the starcode executable
    """
    ssw = ssw_lib.CSsw(ssw_path)

    alphabet, letter_to_int, mat = ssw_build_matrix()

    with pysam.FastaFile(adapter_fasta_filename) as adapter_fasta_file:
        adapter_sequence = adapter_fasta_file.fetch(reference='adapter_sequence')
    with pysam.FastaFile(tso_fasta_filename) as tso_fasta_file:
        tso_sequence = tso_fasta_file.fetch(reference='adapter_sequence')

    with pysam.AlignmentFile(bam_filename, 'rb', check_sq=False) as bam_file:
        with pysam.AlignmentFile(analysis_name + '.intermediate.bam', 'wb', header=bam_file.header) as intermediate_file:
            if whitelist_10x_filename:
                whitelist_10x = read_barcodes(whitelist_10x_filename)
            else:
                whitelist_10x = None

            if whitelist_illumina_filename:
                whitelist_illumina = read_barcodes(whitelist_illumina_filename)
            else:
                whitelist_illumina = None

            stats = AnalysisStats()

            observed_barcodes = dict()
            if record_umis:
                observed_barcodes_umis = dict()

            last_timing = time.time()
            reads_seen = 0

            for read in bam_file.fetch(until_eof=True):
                reads_seen += 1

                if reads_seen % 1000 == 0:
                    current_time = time.time()
                    print('Processed reads: {:,}. Time elapsed: {:.2f}s'.format(reads_seen, current_time - last_timing))
                    last_timing = current_time

                if max_reads is not None and reads_seen > max_reads:
                    break

                stats.reads_seen += 1

                sequence, adapter_alignment_end = align(read, stats, ssw, alphabet, letter_to_int, mat, read_end_length, adapter_sequence, tso_sequence)

                if sequence is None or adapter_alignment_end is None:
                    read.set_tag(ADAPTER_TAG, ".", value_type='Z')
                    read.set_tag(RAW_BARCODE_TAG, ".", value_type='Z')
                    read.set_tag(UMI_TAG, ".", value_type='Z')
                    intermediate_file.write(read)

                    continue

                # Barcode
                observed_barcode, observed_barcode_position, observed_barcode_in_10x, observed_barcode_in_illumina = process_barcode(sequence, adapter_alignment_end, stats, whitelist_10x, whitelist_illumina)
                if observed_barcode is None:
                    read.set_tag(ADAPTER_TAG, adapter_sequence, value_type='Z')
                    read.set_tag(RAW_BARCODE_TAG, ".", value_type='Z')
                    read.set_tag(UMI_TAG, ".", value_type='Z')
                    intermediate_file.write(read)

                    continue

                observed_barcodes[observed_barcode] = observed_barcodes.get(observed_barcode, 0) + 1

                # UMI
                observed_umi, observed_umi_position = process_umi(sequence, observed_barcode_position, stats)
                if observed_umi is None:
                    read.set_tag(ADAPTER_TAG, adapter_sequence, value_type='Z')
                    read.set_tag(RAW_BARCODE_TAG, observed_barcode, value_type='Z')
                    read.set_tag(UMI_TAG, ".", value_type='Z')
                    intermediate_file.write(read)

                    continue

                # Poly-T
                poly_t_found = process_poly_t(sequence, observed_umi_position, stats)
                if not poly_t_found:
                    pass

                observed_barcodes_umis[(observed_barcode, observed_umi)] = observed_barcodes_umis.get((observed_barcode, observed_umi), 0) + 1

                read.set_tag(ADAPTER_TAG, adapter_sequence, value_type='Z')
                read.set_tag(RAW_BARCODE_TAG, observed_barcode, value_type='Z')
                read.set_tag(UMI_TAG, observed_umi, value_type='Z')
                intermediate_file.write(read)

    print('Performing barcode corrections...')

    correct_barcodes(observed_barcodes, analysis_name, stats, whitelist_10x, whitelist_illumina, starcode_path)

    if analysis_name:
        if record_umis:
            with open('{}_barcode_stats.tsv'.format(analysis_name), 'w') as barcode_stats:
                barcode_stats.write(
                    'barcode\tumi\toccurrence\n')
                for ((barcode, umi), occurrence) in observed_barcodes_umis.items():
                    barcode_stats.write('{}\t{}\t{}\n'.format(barcode, umi, int(occurrence)))

        with open('{}_stats.tsv'.format(analysis_name), 'w') as stats_file:
            stats_file.write(AnalysisStats.print_header)
            stats_file.write(stats.print_string())


def correct_barcodes(observed_barcodes, analysis_name, stats, whitelist_10x, whitelist_illumina, starcode_path):
    """
    Performs correction of the observed barcodes and annotates the intermediate file accordingly
    :param observed_barcodes: A dict with the observed barcode sequences as keys and the number of observations as values
    :param analysis_name: Prefix for storing the analysis stats files
    :param stats: The AnalysisStats object
    :param whitelist_10x: The 10x whitelist. Can be None or empty, unless whitelist_illumina_filename is None or empty
    :param whitelist_illumina: The Illumina whitelist. Can be None or empty
    :param starcode_path: Path to the starcode executable
    """
    #subprocess.check_call(["/opt/conda/envs/10x_tool/bin/samtools", "index", analysis_name + '.intermediate.bam'])

    correction_dict = perform_barcode_correction_starcode(observed_barcodes, analysis_name, starcode_path, stats)
    correction_dict['.'] = '.'
    print(len(correction_dict))

    last_timing = time.time()
    reads_in_intermediate_file = 0

    with pysam.AlignmentFile(analysis_name + '.intermediate.bam', 'rb', check_sq=False) as bam_file:
        with pysam.AlignmentFile(analysis_name + '.bam', 'wb', check_sq=False, header=bam_file.header) as output_file:
            for read in bam_file.fetch(until_eof=True):
                if reads_in_intermediate_file % 1000 == 0:
                    current_time = time.time()
                    print('Processed reads: {:,}. Time elapsed: {:.2f}s'.format(reads_in_intermediate_file, current_time - last_timing))
                    last_timing = current_time

                read.set_tag(BARCODE_TAG, read.get_tag(RAW_BARCODE_TAG), value_type='Z')

                reads_in_intermediate_file += 1
                raw_barcode = read.get_tag(RAW_BARCODE_TAG)
                corrected_barcode = correction_dict[raw_barcode]
                if raw_barcode == corrected_barcode:
                    stats.barcode_not_corrected += 1
                else:
                    stats.barcode_corrected += 1

                if whitelist_illumina is not None:
                    if corrected_barcode in whitelist_illumina:
                        read.set_tag(BARCODE_TAG, corrected_barcode, value_type='Z')
                    else:
                        stats.barcode_filtered += 1
                else:
                    if whitelist_10x is not None:
                        if corrected_barcode in whitelist_10x:
                            read.set_tag(BARCODE_TAG, corrected_barcode, value_type='Z')
                        else:
                            stats.barcode_filtered += 1
                    else:
                        read.set_tag(BARCODE_TAG, corrected_barcode, value_type='Z')

                if whitelist_10x is not None:
                    if raw_barcode not in whitelist_10x:
                        if corrected_barcode in whitelist_10x:
                            stats.corrected_from_no_list_to_10x += 1

                output_file.write(read)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads an input BAM file and tries to extract adapter and barcode sequences. When found, annotated reads are written to an output file, if an output filename is provided',
        epilog='Two barcode whitelists can be provided: a 10x whitelist and an Illumina whitelist. If both whitelists are provided, reads will only be annotated if the barcode matches the Illumina list. If only a 10x whitelist is provided, reads will be annotated if the barcode matches the 10x list. In any case, reads will be annotated with the "raw barcode" tag (CR). If no whitelist is provided, all reads will be annotated with the barcode tage (CB) after correction.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-b', '--bam', help='BAM filename', required=True)
    requiredNamed.add_argument('-a', '--adapter', help='Adapter FASTA filename. BWA index must be present.', required=True)
    requiredNamed.add_argument('-r', '--reverse-adapter', help='Reverse adapter FASTA filename. BWA index must be present.', required=True)
    requiredNamed.add_argument('-n', '--name', help='Analysis name (output prefix)', required=True)
    parser.add_argument('--whitelist-10x', help='10x whitelist filename. This may be GZIP compressed (has to have extension .gz in that case)')
    parser.add_argument('--whitelist-illumina', help='Illumina whitelist filename. This may be GZIP compressed (has to have extension .gz in that case)')
    parser.add_argument('--max-reads', help='Number of reads after which the processing should be terminated', type=int, default=None)
    parser.add_argument('--contig', help='Perform analysis only on this contig', default=None)
    parser.add_argument('--read-end-length', help='Interval from both ends of the read in which to search for the adapter sequence', type=int, default=80)
    parser.add_argument('--record-umis', action='store_true', help='If enabled, all barcodes and UMIs will be written to file. This increases memory usage.')
    parser.add_argument('--ssw-path', help='Path to the Striped Smith-Waterman library', type=str, default='/lrma/ssw')
    parser.add_argument('--starcode-path', help='Path to the starcode executable', type=str, default='/lrma/starcode-master/starcode')

    args = parser.parse_args()

    if args.whitelist_illumina and not args.whitelist_10x:
        print('Illumina whitelist provided but no 10x whitelist provided.')
        exit(1)

    main(args.bam, args.name, args.adapter, args.reverse_adapter, args.whitelist_10x, args.whitelist_illumina, args.max_reads, args.contig, args.read_end_length, args.record_umis, args.ssw_path, args.starcode_path)
