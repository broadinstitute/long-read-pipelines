#!/usr/bin/env python3.8

import time
import argparse
import math
import logging
import os
import sys
import inspect
import subprocess
import tempfile
import re
import uuid

from collections import OrderedDict
from collections import namedtuple

from enum import Enum

import numpy as np

import pysam

################################################################################


USE_BWA_MEM_2 = False

RC_READ_NAME_IDENTIFIER = "_RC"
VALID_DNA_SEQUENCE_PATTERN = re.compile(r"^[ATGCNatgcn]+$")

# Default max PL:
MAX_ALIGNMENT_PL = 60

# Minimum reported PL quality value (will override lower values):
MIN_ALIGNMENT_PL = 0

# Min PL for an alignment to be kept as "good":
# This should be equal to the base read quality.
MIN_GOOD_ALIGNMENT_PL = 7.0

# Min number of bases required to be in an alignment for it to be kept in consideration.
# TODO: Make the default based on the number of bases in each known segment.
MIN_ALIGNMENT_LENGTH = 4

# Named tuple to store alignment information:
DetailedAlignmentInfo = namedtuple(
    "DetailedAlignmentInfo", ["start_pos", "end_pos", "template_length", "cigar", "qual_pl"]
)

ReadNameAndSeq = namedtuple(
    "ReadNameAndSeq", ["name", "seq", "raw_obj"]
)


class ProcessedAlignmentResult(namedtuple(
    "ProcessedAlignmentResult", ["seq_name", "alignment_string", "target_start_index", "target_end_index",
                                 "read_start_pos", "read_end_pos", "template_length", "cigar", "overall_quality"]
)):
    __slots__ = ()

    def __str__(self):

        return (
            f"ProcessedAlignmentResult({self.seq_name}, {self.alignment_string}, {self.target_start_index}, "
            f"{self.target_end_index}, {self.read_start_pos}, {self.read_end_pos}, "
            f"{self.template_length}, {cigar_tuple_to_string(self.cigar)}, {self.overall_quality})"
        )


################################################################################

LOGGER = logging.getLogger("extract_bounded_read_sections")


def configure_logging(args):
    """Set up logging for the module"""

    format_string = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'

    # Set logging level:
    log_level = logging.INFO
    if args.quiet:
        log_level = logging.CRITICAL
    elif args.verbose:
        log_level = logging.DEBUG
    elif args.veryverbose:
        log_level = logging.NOTSET

    logging.basicConfig(level=log_level, format=format_string)


def print_logo():
    """Print the logo to the log."""

    LOGGER.info("====================================================================")
    LOGGER.info("    ____           _                              _")
    LOGGER.info("   / ___|__ _ _ __| |_ ___   __ _ _ __ __ _ _ __ | |__   ___ _ __")
    LOGGER.info("  | |   / _` | '__| __/ _ \\ / _` | '__/ _` | '_ \\| '_ \\ / _ \\ '__|")
    LOGGER.info("  | |__| (_| | |  | || (_) | (_| | | | (_| | |_) | | | |  __/ |")
    LOGGER.info("   \\____\\__,_|_|   \\__\\___/ \\__, |_|  \\__,_| .__/|_| |_|\\___|_|")
    LOGGER.info("                            |___/          |_|")
    LOGGER.info(" ___________________________________________________________________")
    LOGGER.info("")
    LOGGER.info("                   Split Sequence on Delimiters")
    LOGGER.info(" ___________________________________________________________________")
    LOGGER.info("/                   _____.----------._                              \\")
    LOGGER.info("|  _         __,---'_       \"         `.     _,--.       __,---.    |")
    LOGGER.info("| | \\___    /      ( )   \"          \"   `-.,' (') \\     (  / _\\ \\   |")
    LOGGER.info("|  \\    \\__/_   __(( _)_      (    \"   \"     (_\\_) \\___  `-.___,'   |")
    LOGGER.info("|   \\     (  )_(__)_|( ))  \"   ))          \"   |    \"  \\        _   |")
    LOGGER.info("|    \\__ (( _( (  ))  ) _)    ((     \\\\//    \" |   \"    \\_____,' |  |")
    LOGGER.info("|       \\  ( ))(_)(_)_)|  \"    ))    //\\\\ \" __,---._  \"  \"   \"  /   |")
    LOGGER.info("|        |(_ _)| | |   |   \"  (   \"      ,-'        `-.   ___  /    |")
    LOGGER.info("|        |  |  |   |   _,--- ,--. _  \"  (              ) /___\\ \\    |")
    LOGGER.info("|       /   |      _,----._,'`--'\\.`-._  `._  _ __ _,-'  |H__|  \\   |")
    LOGGER.info("|      / \"     _,-' / `\\ ,' / _'  \\`.---.._          __        \" \\  |")
    LOGGER.info("|     / /   .-' , / ' _,'_  -  _ '- _`._ `.`-._    _/- `--.   \" \" \\ |")
    LOGGER.info("|    / / _-- `---, .-' __   --  _,---.  `-._   _,-'- / ` \\ \\_   \" | |")
    LOGGER.info("|   | | -- _    / /  `-_- _  _,' '  \\ \\_`-._,-'  / --   \\  - \\_   / |")
    LOGGER.info("|   | \\ -      /  | \"     ,-'_ /-  `_ ._`._`-...._____...._,--'  /  |")
    LOGGER.info("|   \\  \\_ /   /  /    ___  `---  ---  - - ' ,--.     ___        |   |")
    LOGGER.info("|    \\      ,'  |  \" (o o)   \"         \" \" |    \\_,-'   `.     ,'   |")
    LOGGER.info("|     |__,-'     \\    \\\"/      \"  \"   \"    /       O     `-.__/     |")
    LOGGER.info("|                 `.______________________/        |                |")
    LOGGER.info("\\___________________________________________________________________/")
    LOGGER.info("")
    LOGGER.info("")
    if USE_BWA_MEM_2:
        LOGGER.info("      powered by BWA MEM 2")
    else:
        LOGGER.info("      powered by BWA MEM")
    LOGGER.info("====================================================================")
    LOGGER.info("")


################################################################################


class IdentityMap:
    """A map to that returns exactly what is given."""
    def __init__(self):
        pass

    def __getitem__(self, key):
        return key


class ReadFile:
    """
    A class to abstract away the differences between reading from bam/sam files and fasta/fastq files with pysam.
    """
    def __init__(self, file_name):
        self.file_name = file_name

        self._file_object = None
        self._is_alignment_file = file_name.lower().endswith(".sam") or file_name.lower().endswith(".bam")

    def __enter__(self):

        if self._is_alignment_file:
            # Determine read flags (is it a bam or a sam file?)
            file_flags = ''
            if self.file_name.endswith('.bam'):
                file_flags = 'b'

            self._file_object = pysam.AlignmentFile(self.file_name, 'r' + file_flags, check_sq=False)
        else:
            self._file_object = pysam.FastxFile(self.file_name)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file_object.close()

    def is_sam(self):
        """
        :return: True iff this file is a sam/bam file.
        """
        return self._is_alignment_file

    def get_header(self):
        """:return: The header of the underlying alignment file or None."""
        return self._file_object.header if self._is_alignment_file else None

    def get_reads(self):
        """
        Generator function to yield a Sequence object for every read in the
        supporting read file corresponding to self._file_name
        """

        if self._is_alignment_file:
            for read in self._file_object.fetch(until_eof=True):
                yield ReadNameAndSeq(read.query_name, read.query_sequence, read)
        else:
            for read in self._file_object:
                yield ReadNameAndSeq(read.name, read.sequence, read)


# IUPAC RC's from: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
# and https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
RC_BASE_MAP = {"N": "N", "A": "T", "T": "A", "G": "C", "C": "G", "Y": "R", "R": "Y", "S": "S", "W": "W", "K": "M",
               "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "n": "n", "a": "t", "t": "a", "g": "c", "c": "g",
               "y": "r", "r": "y", "s": "s", "w": "w", "k": "m", "m": "k", "b": "v", "v": "b", "d": "h", "h": "d"}


def reverse_complement(base_string):
    """
    Reverse complements the given base_string.
    :param base_string: String of bases to be reverse-complemented.
    :return: The reverse complement of the given base string.
    """

    return ''.join(map(lambda b: RC_BASE_MAP[b], base_string[::-1]))


CIGAR_ELEMENT_STRING_MAP = {
    pysam.CDEL: "D",
    pysam.CMATCH: "M",
    pysam.CBACK: "B",
    pysam.CINS: "I",
    pysam.CPAD: "P",
    pysam.CEQUAL: "=",
    pysam.CSOFT_CLIP: "S",
    pysam.CHARD_CLIP: "H",
    pysam.CREF_SKIP: "N",
    pysam.CDIFF: "X"
}


def cigar_tuple_to_string(cigar_tuples):
    """
    Converts a given list of cigar tuples to a string.
    :param cigar_tuples: list of tuples, each containing a cigar element and count.
    :return: A string representing the given cigar tuple list.
    """

    cigar_string_elements = []
    for element, count in cigar_tuples:
        cigar_string_elements.append(f"{count}{CIGAR_ELEMENT_STRING_MAP[element]}")

    return "".join(cigar_string_elements)


def get_pg_sam_header_entry():
    """Returns a dictionary containing the PG entry for the invocation of this file.
    This returned value is suitable to adding directly to an output SAM header dictionary."""

    name = sys.argv[0][sys.argv[0].rfind(os.path.sep):]

    return {"PG": {
        "ID": name,
        "PN": name,
        "VN": "UNKNOWN",
        "DS": "Splits input reads by given delimiters into smaller read sequences.",
        "CL": " ".join(sys.argv),
    }}


def _log_var(var):
    """
    Logs the given variable's name and contents at the DEBUG log level.

    Supports logging even if var is an expression in the invocation.
    :param var: The variable whose contents are to be logged.
    :return: None
    """
    # prev_frame = inspect.currentframe().f_back

    # Get the name of this function:
    fn_name = inspect.currentframe().f_code.co_name

    # Get the source code line of where this function was called:
    call_line = inspect.stack()[1][4][0].strip()

    # Make sure we're in the right place:
    assert call_line.startswith(f"{fn_name}(")

    # Pull out the variable name from the function call:
    var_name = call_line[len(f"{fn_name}("):][:-1].strip()

    LOGGER.debug("%s = %s", var_name, var)


################################################################################


def dump_seq_map(seq_map, name='Reads'):
    """Dumps the given sequence map to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug('%s:', name)
        for e in seq_map.items():
            LOGGER.debug('    %s -> %s', e[0], e[1])


def assert_unique_delimiter_seqs(delim_name_map):
    """Asserts that the sequences of the given delimiters are unique."""
    seq_set = set()
    for k,v in delim_name_map.items():
        if v in seq_set:
            LOGGER.error(f"ERROR: Delimiter sequence is not unique for seq {k}: {v}")
            LOGGER.error("ABORTING.")
            sys.exit(1)
        else:
            seq_set.add(v)


def ingest_fastx_file(file_path, allow_duplicate_names=False):
    """Ingest the contents of a FASTA/FASTQ file and return two dictionaries
    
        1: Mapping from read name to read sequence
        2: Mapping from template name to read name

    The `template name` is simply the word 'template' with the
    order in which a given read occurs in the file
    (e.g. 'template0' or 'template10').
    """

    t_num = 0
    _read_to_sequence_dict = OrderedDict()
    _template_to_read_name_dict = dict()

    _name_set = set()

    with pysam.FastxFile(file_path) as file_handle:
        for entry in file_handle:

            if not allow_duplicate_names and entry.name in _name_set:
                LOGGER.error(f"ERROR: Delimiter name is not unique: {entry.name}")
                LOGGER.error("ABORTING.")
                sys.exit(1)

            _read_to_sequence_dict[entry.name] = entry.sequence
            _template_to_read_name_dict[f"template{t_num}"] = entry.name

            _name_set.add(entry.name)

            t_num += 1

    return _read_to_sequence_dict, _template_to_read_name_dict


def get_qual_pl(num_errors, seq_length):
    """
    Computes the Phred-Scaled quality score of an alignment.
    :param num_errors: Number of errors in the alignment.
    :param seq_length: Length of the alignment.
    :return: The PL quality.
    """

    if num_errors == 0:
        return MAX_ALIGNMENT_PL
    else:
        q = -10 * math.log10(num_errors / seq_length)
        if q < MIN_ALIGNMENT_PL:
            return MIN_ALIGNMENT_PL
        else:
            return q


def create_alignment_with_bwa_mem(read_data,
                                  delimiter_names_to_seq_dict,
                                  minqual,
                                  minbases,
                                  threads=1,
                                  log_spacing="    "):
    """
    Align the given delimiters to the given read_data using BWA MEM 2.

    :param read_data: ReadNameAndSeq object against which to align all sequences in target_sequences
    :param delimiter_names_to_seq_dict: dictionary of delimiter names to sequences
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :param threads: number of threads to use when aligning.
    :param log_spacing: Spacing to precede any log statements.
    :return: A list of ProcessedAlignmentResult objects.
    """
    out_file_name = "tmp.sam"

    bwa_exe_path = "/bwa-mem2-2.0pre2_x64-linux/bwa-mem2" if USE_BWA_MEM_2 else "/bwa/bwa"

    # Write sequences to tmp fasta file:
    _, ref_file = tempfile.mkstemp()
    _, delim_file = tempfile.mkstemp()
    try:
        LOGGER.debug("%sCreating tmp \"reference\" fasta file: %s", log_spacing, ref_file)
        with open(ref_file, "w", encoding="ascii") as tmp:
            tmp.write(f">{read_data.name}\n")
            tmp.write(f"{read_data.seq}\n")

        bwa_index_args = [bwa_exe_path, "index", ref_file]
        LOGGER.debug(
            "%sRunning BWA Index on tmp file (%s): %s", log_spacing, ref_file, " ".join(bwa_index_args)
        )
        _ = subprocess.run(bwa_index_args, capture_output=True, check=True)

        LOGGER.debug("%sCreating tmp delimiter fasta file: %s", log_spacing, delim_file)
        with open(delim_file, "w", encoding="ascii") as tmp:
            for delim_name, delim_seq in delimiter_names_to_seq_dict.items():
                tmp.write(f">{delim_name}\n{delim_seq}\n")

        LOGGER.debug("Contents of tmp \"reference\" fasta file:")
        with open(ref_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        LOGGER.debug("Contents of known segment \"read\" fasta file:")
        with open(delim_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        bwa_mem_args = [bwa_exe_path, "mem",
                        "-a",  # Output all found alignments for single-end or unpaired paired-end reads.
                        # These alignments will be flagged as secondary alignments.
                        "-S",  # skip mate rescue
                        "-P",  # skip pairing; mate rescue performed unless -S also in use
                        "-k8",  # minimum seed length
                        "-A", "1",  # Matching score.
                        "-B", "4",  # Mismatch penalty.
                        #  The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}.
                        "-O", "6,6",  # Gap open penalty.
                        "-E", "1,1",  # gap extension penalty; a gap of size k cost '{-O} + {-E}*k'
                        "-L", "5,5",  # penalty for 5'- and 3'-end clipping
                        "-U", "17",  # penalty for an unpaired read pair
                        "-T", "30",  # minimum score to output
                        "-c", "1000",  # skip seeds with more than INT occurrences
                        "-t", str(threads),
                        "-o", out_file_name,
                        ref_file, delim_file]
        LOGGER.debug(
            "%sRunning BWA Mem on tmp file (%s): %s", log_spacing, ref_file, " ".join(bwa_mem_args)
        )
        completed_process = subprocess.run(bwa_mem_args,
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)

        if LOGGER.isEnabledFor(logging.DEBUG):
            LOGGER.debug("BWA Mem output:")
            for l in completed_process.stdout.split("\n"):
                LOGGER.debug(l)

            with open(out_file_name, 'rb') as f:
                LOGGER.debug("=" * 80)
                LOGGER.debug("Raw BWA MEM Alignment:")
                for line in f.readlines():
                    LOGGER.debug("%s", line.decode("ascii").rstrip())
                LOGGER.debug("=" * 80)

        return get_processed_results_from_bwa_mem_file(out_file_name, minqual, minbases)

    except subprocess.CalledProcessError as e:
        LOGGER.error("Could not align with BWA Mem!")
        LOGGER.error("Stdout: %s", e.stdout.decode("utf-8"))
        LOGGER.error("Stderr: %s", e.stderr.decode("utf-8"))
        raise e

    finally:
        os.remove(ref_file)
        os.remove(delim_file)
        try:
            os.remove(out_file_name)
            pass
        except FileNotFoundError:
            # If the alignment failed, we won't necessarily have an output file:
            pass


def get_processed_results_from_bwa_mem_file(file_path, minqual, minbases):
    """
    Ingests the given file and creates raw alignments from it.
    :param file_path: Path to the sam file of a BWA MEM run.
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :return: A list of ProcessedAlignmentResult objects.
    """

    processed_results = []

    # Only primary hits will have sequence information with them so we have to
    # Read it off first.  There shouldn't be very many reads, so this is a little slow, but should be OK.
    read_seqs = dict()
    with pysam.AlignmentFile(file_path, 'r', check_sq=False) as f:
        for read in f.fetch(until_eof=True):
            if read.query_sequence:
                read_seqs[read.query_name] = read.query_sequence

    with pysam.AlignmentFile(file_path, 'r', check_sq=False) as f:
        for read in f.fetch(until_eof=True):

            if read.is_unmapped:
                continue

            seq_name = read.query_name
            if read.is_reverse:
                seq_name = seq_name + RC_READ_NAME_IDENTIFIER

            bases = read_seqs[read.query_name]

            template_length = read.infer_read_length()
            qual_pl = get_qual_pl(read.get_tag("NM"), template_length)

            # Get the leading and trailing clips so we can remove them from the aligned string:
            leading_clips = 0
            for e in read.cigartuples:
                if e[0] == pysam.CSOFT_CLIP or e[0] == pysam.CHARD_CLIP:
                    leading_clips += e[1]
                else:
                    break

            trailing_clips = 0
            for e in read.cigartuples[::-1]:
                if e[0] == pysam.CSOFT_CLIP or e[0] == pysam.CHARD_CLIP:
                    trailing_clips += e[1]
                else:
                    break

            # Note - must adjust ref start/end pos to align properly with conventions from other aligners
            #        (other aligner conventions: 1-based coordinates, inclusive end positions)
            p = ProcessedAlignmentResult(
                seq_name, bases, leading_clips, len(bases)-trailing_clips-1,
                int(read.reference_start), int(read.reference_end-1),
                template_length, tuple(read.cigartuples), qual_pl
            )

            # Check against thresholds to make sure we should report the alignment:
            if qual_pl < minqual or template_length < minbases:
                if qual_pl < minqual and template_length < minbases:
                    reason_string = f"qual too low ({qual_pl} < {minqual}) " \
                                    f"AND aligment too short ({template_length} < {minbases})"
                elif template_length < minbases:
                    reason_string = f"aligment too short ({template_length} < {minbases})"
                else:
                    reason_string = f"qual too low ({qual_pl} < {minqual})"

                LOGGER.debug("Target does not pass threshold: %s: %s (%s)", seq_name, reason_string, p)
                continue

            # Check for secondary alignment flags:
            if read.is_secondary:
                LOGGER.debug("Skipping - Target is a secondary alignment: %s", p)
                continue
            elif read.is_supplementary:
                LOGGER.debug("Skipping - Target is a supplementary alignment: %s", p)
                continue

            # Read passes all filters.  We can add it to our list:
            processed_results.append(p)

    # Sort by the order in which they appear in the read:
    processed_results.sort(key=lambda x: x.read_start_pos)

    for r in processed_results:
        LOGGER.debug("    %s", r)

    return processed_results


def filter_alignment_results_by_position(segment_alignment_results):
    """Filter the given alignment results to contain only one delimiter at each read start position."""
    filtered_results = []

    LOGGER.info(f"Applying filtering {len(segment_alignment_results)} results...")

    # This isn't the best way to do this - really it should happen down in the alignment stage
    # so we don't have to iterate so many times.
    segment_dict = dict()
    for s in segment_alignment_results:
        if s.read_start_pos in segment_dict:
            segment_dict[s.read_start_pos].append(s)
        else:
            segment_dict[s.read_start_pos] = [s]

    # Now that we have our map we can perform the filtering:
    num_filtered = 0
    for k,v in segment_dict.items():
        # Most of the time we should have only one alignment per segment:
        if len(v) == 1:
            LOGGER.debug(f"Single sequence detected for position: {k}")
            filtered_results.append(v[0])
        else:
            # Since we have more than one alignment,
            # we need to get the alignment with the best score:
            v.sort(key=lambda x: x.overall_quality)
            filtered_results.append(v[0])

            if LOGGER.isEnabledFor(logging.DEBUG):
                LOGGER.debug(f"Filtering {len(v) - 1} results for position: {k}: %s",
                             ",".join([f"{s.seq_name}" for s in v[1:]]))

            num_filtered += len(v) - 1

    LOGGER.info(f"Filtered {num_filtered} results.")
    LOGGER.debug(f"Returning {len(filtered_results)} results.")

    return filtered_results


def align_delimiters(read_data, delimiter_names_to_seq_dict, minqual, minbases, threads=1, log_spacing="    "):
    """
    Align the given delimiters to the given read sequence.

    :param read_data: ReadNameAndSeq object against which to align all sequences in target_sequences
    :param delimiter_names_to_seq_dict: dictionary of delimiter names to sequences
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :param threads: Number of threads with which to run Tesserae alignment.
    :param log_spacing: Spacing to precede any log statements.
    :return: A list of ProcessedAlignmentResult objects.
    """

    start_time = time.time()
    processed_results = create_alignment_with_bwa_mem(read_data, delimiter_names_to_seq_dict, minqual, minbases,
                                                      threads=threads)
    end_time = time.time()
    LOGGER.info("%sCreated %d Alignments.  Alignment took %fs", log_spacing, len(processed_results),
                end_time - start_time)

    return processed_results


def write_sub_sequences(read_data, aligned_delimiters, out_bam_file, out_tsv_file):
    """
    Write out the sections of the given read_data as bounded by aligned_delimiters to the given out_file in
    BAM format.
    Also write out the segment splitting information in TSV format in the given TSV output file.
    :param read_data: A ReadNameAndSeq object representing the parent read to the given aligned_delimiters.
    :param aligned_delimiters: A list of ProcessedAlignmentResult objects.
    :param out_bam_file: An open pysam.AlignmentFile object to which to write results.
    :param out_tsv_file: An open file object to which to write alignments in TSV format.
    :return: None
    """

    cur_read_base_index = 0
    prev_delim_name = "START"

    out_tsv_file.write(f"{read_data.name}\t{len(read_data.seq)}\t")

    for i, delimiter_alignment in enumerate(aligned_delimiters):

        # Do some math here to account for how we create the tuple list:
        # And adjust for reverse complements:
        start_coord = cur_read_base_index
        end_coord = delimiter_alignment.read_start_pos
        delim_name = delimiter_alignment.seq_name

        # Write out bam:
        a = pysam.AlignedSegment()
        a.query_name = f"{read_data.name}_{start_coord+1}-{end_coord}_{prev_delim_name}-{delim_name}"
        a.query_sequence = f"{read_data.seq[start_coord:end_coord]}"
        a.query_qualities = read_data.raw_obj.query_alignment_qualities[start_coord:end_coord]
        a.tags = read_data.raw_obj.get_tags()
        a.flag = 4  # unmapped flag
        a.mapping_quality = 255
        out_bam_file.write(a)

        # Write the TSV output:
        out_tsv_file.write(f"{delimiter_alignment.seq_name}"
                           f"[{delimiter_alignment.read_start_pos}-{delimiter_alignment.read_end_pos}]"
                           f"c{cigar_tuple_to_string(delimiter_alignment.cigar)}"
                           f"@{delimiter_alignment.overall_quality:2.3f}")
        if i != len(aligned_delimiters) - 1:
            out_tsv_file.write(",")

        cur_read_base_index = end_coord
        prev_delim_name = delim_name

    # Now we have to write out the last segment:
    start_coord = cur_read_base_index
    end_coord = len(read_data.seq)
    delim_name = "END" 
    
    a = pysam.AlignedSegment()
    a.query_name = f"{read_data.name}_{start_coord+1}-{end_coord}_{prev_delim_name}-{delim_name}"
    a.query_sequence = f"{read_data.seq[start_coord:end_coord]}"
    a.query_qualities = read_data.raw_obj.query_alignment_qualities[start_coord:end_coord]
    a.tags = read_data.raw_obj.get_tags()
    a.flag = 4  # unmapped flag
    a.mapping_quality = 255
    out_bam_file.write(a)

    out_tsv_file.write('\n')


def split_sequences(args):
    """Main CLI call for the Extract Bounded Read Sections tool."""

    # Set up multi-threadding if the system will support it:
    LOGGER.info("Python version: %s", sys.version.replace('\n', ''))
    min_multithread_version = (3, 8)
    if sys.version_info >= (3, 8):
        num_threads = os.cpu_count() - 1
        num_threads = 1 if num_threads < 1 else num_threads
    else:
        LOGGER.warning(
            "Python version to early for multithreading (%d.%d.%d<%d.%d).",
            sys.version_info[0],
            sys.version_info[1],
            sys.version_info[2],
            min_multithread_version[0],
            min_multithread_version[1],
        )
        num_threads = 1
    LOGGER.info("Setting thread count to: %d", num_threads)

    if args.max_read_length:
        LOGGER.info("Filtering out reads of length > %d to file: %s", args.max_read_length, args.rejected_outfile)

    LOGGER.info("Writing output to %s", args.outfile)
    if os.path.exists(args.outfile):
        LOGGER.warning("Outfile already exists.  Will overwrite: %s", args.outfile)

    LOGGER.info("Ingesting delimiters from %s ...", args.delimiters)
    delimiter_names_to_seq_dict, _ = ingest_fastx_file(args.delimiters)
    assert_unique_delimiter_seqs(delimiter_names_to_seq_dict)
    LOGGER.info("Ingested %d delimiters.", len(delimiter_names_to_seq_dict))
    dump_seq_map(delimiter_names_to_seq_dict, "Delimiters")

    # A spacing variable for nice looking logs:
    spacing_one = " " * 4

    # Open all our files here so they'll be automatically closed:
    with ReadFile(args.reads) as reads_file:

        out_bam_header_dict = reads_file.get_header().to_dict()
        if out_bam_header_dict is None:
            out_bam_header_dict = {'HD': {'VN': '1.0', 'SO': "unknown", 'pb': '>=3.01'}}

        with pysam.AlignmentFile(args.rejected_outfile, 'wb', header=out_bam_header_dict) as rejected_out_file, \
                pysam.AlignmentFile(args.outfile, 'wb', header=out_bam_header_dict) as out_bam, \
                open(args.annotated_read_sections_file, 'w') as read_segments_file:

            num_delimiters_detected = 0
            num_reads_with_delimiters = 0
            num_forward_subsequences_extracted = 0
            num_rc_subsequences_extracted = 0
            num_rejected = 0

            LOGGER.debug("Writing TSV split file header...")
            read_segments_file.write("\t".join(["READ_NAME", "READ_LENGTH", "SEGMENTS"]))
            read_segments_file.write("\n")

            LOGGER.info("Processing reads...")

            num_reads = 0
            for read_num, read_data in enumerate(reads_file.get_reads()):
                num_reads += 1

                LOGGER.debug(
                    "%sProcessing read %d: %s (len=%d)",
                    spacing_one,
                    read_num,
                    read_data.name,
                    len(read_data.seq)
                )

                if args.max_read_length and len(read_data.seq) > args.max_read_length:
                    LOGGER.warning("Ignoring read %d - %s: Length too long: %d > %d", read_num, read_data.name,
                                   len(read_data.seq), args.max_read_length)
                    rejected_out_file.write(read_data.raw_obj)
                    num_rejected += 1
                    continue

                LOGGER.info("%sAlignment of delimiters ...", spacing_one)
                segment_alignment_results = align_delimiters(
                    read_data, delimiter_names_to_seq_dict, args.minqual, args.minbases, threads=num_threads
                )

                if len(segment_alignment_results) == 0:
                    LOGGER.info("%sThis read contains no delimiters.", spacing_one)
                else:
                    LOGGER.info("%sDelimiters occur %d time(s).", spacing_one, len(segment_alignment_results))

                    # Filter the delimiters down to those we'll keep:
                    filtered_alignment_results = filter_alignment_results_by_position(segment_alignment_results)

                    # Write out our new subsequences to the output fasta file:
                    write_sub_sequences(read_data, filtered_alignment_results, out_bam, read_segments_file)

                    # Track subsequence statistics:
                    num_forward_subsequences_extracted += sum(
                        not b.seq_name.endswith(RC_READ_NAME_IDENTIFIER) for b in filtered_alignment_results
                    )
                    num_rc_subsequences_extracted += sum(
                        b.seq_name.endswith(RC_READ_NAME_IDENTIFIER) for b in filtered_alignment_results
                    )
                    num_delimiters_detected_this_read = len(filtered_alignment_results) + 1
                    num_delimiters_detected += num_delimiters_detected_this_read
                    num_reads_with_delimiters += 1

            LOGGER.info("Processed %d reads.", num_reads)
            if args.max_read_length:
                LOGGER.info("Rejected %d reads (reads longer than %d bases).", num_rejected, args.max_read_length)
            LOGGER.info("# Reads containing delimiters: %d", num_reads_with_delimiters)
            LOGGER.info("# forward direction delimiters detected: %d", num_forward_subsequences_extracted)
            LOGGER.info("# reverse-complemented direction delimiters detected: %d",
                        num_rc_subsequences_extracted)
            LOGGER.info("Total # sub-sequences extracted: %d", num_delimiters_detected)
            LOGGER.info("Overall data multiplication rate: %2.2f", num_delimiters_detected / num_reads)


################################################################################

def main(raw_args):

    base_outfile_name = "split_sequences"

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        description="Ingests two files: "
                    "one read file(FASTA/FASTQ/SAM/BAM) containing reads to be mapped, "
                    "one FASTA file containing known possible sequences that can occur in the read"
                    "Under the hood uses BWA-MEM 2",
        usage="Naively split every sequence in the given reads file by the delimiters in the given delimters file.",
    )

    align_required_args = parser.add_argument_group("required arguments")
    align_required_args.add_argument(
        "-r", "--reads", help="Reads SAM/BAM/FASTA/FASTQ file.", required=True
    )
    align_required_args.add_argument(
        "-d", "--delimiters", help="Delimiters FASTA/FASTQ file.", required=True
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="Output bam file in which to store alignment results.",
        default=f"{base_outfile_name}.bam",
        required=False,
    )

    parser.add_argument(
        "-m",
        "--minqual",
        help=f"Minimum quality for good alignment (default: {MIN_GOOD_ALIGNMENT_PL})",
        default=MIN_GOOD_ALIGNMENT_PL,
        type=float
    )

    parser.add_argument(
        "-n",
        "--minbases",
        help=f"Minimum number of bases for an alignment to be retained (default: {MIN_ALIGNMENT_LENGTH})",
        default=MIN_ALIGNMENT_LENGTH,
        type=float
    )

    parser.add_argument(
        "--max_read_length",
        help=f"Maximum read length to be processed.  Reads exceeding this length will be written to a rejects file.  "
             f"This option is off by default (no filtering will occur).",
        type=int,
        required=False,
    )

    parser.add_argument(
        "--rejected_outfile",
        help="Output bam file in which to store rejected reads.",
        default=f"{base_outfile_name}.rejected.bam",
        required=False,
    )

    parser.add_argument(
        "--annotated_read_sections_file",
        help="Output TSV in which to put read segmentation information.",
        default=f"{base_outfile_name}.annotations.tsv",
        required=False,
    )

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "-q", "--quiet", help="silence logging except errors", action="store_true"
    )
    verbosity_group.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )
    verbosity_group.add_argument(
        "-vv", "--veryverbose", help="maximal output verbosity", action="store_true"
    )

    # ---------------------------------------

    # Parse args
    args = parser.parse_args()

    configure_logging(args)

    # Print logo:
    print_logo()

    # Log our command-line and log level so we can have it in the log file:
    LOGGER.info("Invoked by: %s", " ".join(raw_args))
    LOGGER.info("Complete runtime configuration settings:")
    for name, val in vars(args).items():
        LOGGER.info("    %s = %s", name, val)
    LOGGER.info("Log level set to: %s", logging.getLevelName(logging.getLogger().level))

    # Call our main method:
    split_sequences(args)

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == '__main__':
    main(sys.argv)
