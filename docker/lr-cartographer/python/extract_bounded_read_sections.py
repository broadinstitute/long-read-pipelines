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

from tesserae import Tesserae
from tesserae.sequence import Sequence
from tesserae.tesserae import TesseraeAlignmentResult
from tesserae.tesserae import DEFAULT_REC

################################################################################

LOGGER = logging.getLogger("extract_bounded_read_sections")

################################################################################

BARCODE_NAME_IDENTIFIER = "BARCODE"
UNKNOWN_NAME_IDENTIFIER = "UNKNOWN"

# Make a random string at the end to attempt to prevent unplanned matches:
RC_READ_NAME_IDENTIFIER_BASE = "_RC"
RC_READ_NAME_IDENTIFIER = RC_READ_NAME_IDENTIFIER_BASE + "_" + str(uuid.uuid1())[:8]


class AlignmentAlgorithm(Enum):
    """
    Specifier for an alignment algorithm.
    """
    TESSERAE = 1
    MOSAIC_ALIGNER = 2
    SMITH_WATERMAN = 3
    NEEDLEMAN_WUNSCH = 4
    BWA_MEM = 5
    BWA_ALN = 6


VALID_DNA_SEQUENCE_PATTERN = re.compile(r"^[ATGCNatgcn]+$")
MOSAIC_ALIGNER_LINE_MATCH_PATTERN = re.compile(r"^\s*(.*?)\s\[\s*(\d+)-\s*(\d+)\]\t+(.*)\s*$")

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

P_REC_KNOWN = DEFAULT_REC
P_REC_UNKNOWN = DEFAULT_REC / 10

# Named tuple to store alignment information:
DetailedAlignmentInfo = namedtuple(
    "DetailedAlignmentInfo", ["start_pos", "end_pos", "template_length", "cigar", "qual_pl"]
)


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

    def get_reads(self):
        """
        Generator function to yield a Sequence object for every read in the
        supporting read file corresponding to self._file_name
        """

        if self._is_alignment_file:
            for read in self._file_object.fetch(until_eof=True):
                yield Sequence(read.query_name, read.query_sequence)
        else:
            for read in self._file_object:
                yield Sequence(read.name, read.sequence)


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
STRING_CIGAR_ELEMENT_MAP = {v: k for (k, v) in CIGAR_ELEMENT_STRING_MAP.items()}


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


MetaSequenceInfo = namedtuple(
    "MetaSequenceInfo", ["name", "raw_read_alignment_string", "alignment_start_index", "alignment_end_index"]
)

################################################################################


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


def print_logo(alignment_type):
    """Print the logo to the log."""

    alignment_algorithm = AlignmentAlgorithm.TESSERAE
    if alignment_type:
        try:
            alignment_algorithm = AlignmentAlgorithm[alignment_type]
        except KeyError:
            LOGGER.error("Error: You must provide a valid alignment algorithm.  Options are: %s",
                         ", ".join([e.name for e in AlignmentAlgorithm]))
            sys.exit(1)

    LOGGER.info("====================================================================")
    LOGGER.info("    ____           _                              _")
    LOGGER.info("   / ___|__ _ _ __| |_ ___   __ _ _ __ __ _ _ __ | |__   ___ _ __")
    LOGGER.info("  | |   / _` | '__| __/ _ \\ / _` | '__/ _` | '_ \\| '_ \\ / _ \\ '__|")
    LOGGER.info("  | |__| (_| | |  | || (_) | (_| | | | (_| | |_) | | | |  __/ |")
    LOGGER.info("   \\____\\__,_|_|   \\__\\___/ \\__, |_|  \\__,_| .__/|_| |_|\\___|_|")
    LOGGER.info("                            |___/          |_|")
    LOGGER.info(" ___________________________________________________________________")
    LOGGER.info("")
    LOGGER.info("                   Extract Bounded Read Sections")
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

    if alignment_algorithm == AlignmentAlgorithm.MOSAIC_ALIGNER:
        LOGGER.info(r"                                   +---------------------------+")
        LOGGER.info(r"                                   |__   \__/  \__   \__/  \__ |")
        LOGGER.info(r"                 powered           |  \__/  \     \__/  \     \| ")
        LOGGER.info(r"                                   |__/     /   __/     /   __/|")
        LOGGER.info(r"                    by             |  \__   \__/  \__   \__/  \| ")
        LOGGER.info(r"                                   |     \__/  \     \__/  \   |")
        LOGGER.info(r"                     MosaicAligner |   __/     /   __/     /   |")
        LOGGER.info(r"                                   |__/  \__   \__/  \__   \__/|")
        LOGGER.info(r"                                   |  \     \__/  \     \__/  \| ")
        LOGGER.info(r"                                   +---------------------------+")
    elif alignment_algorithm == AlignmentAlgorithm.TESSERAE:
        LOGGER.info("                                                     +-----+")
        LOGGER.info("                                                     | \\   | \\")
        LOGGER.info("                                           powered   |  +-----+")
        LOGGER.info("                                                     +--|--+  |")
        LOGGER.info("                                               by     \\ |   \\ |")
        LOGGER.info("                                                        +-----+")
        LOGGER.info("                                                 tesserae")
    else:
        LOGGER.info("")
        LOGGER.info("      powered by %s", alignment_algorithm.name)

    LOGGER.info("====================================================================")
    LOGGER.info("")


def compute_detailed_alignment_info(
    query_alignment_string, target_alignment_string, target_start_index, target_end_index,
):
    """Compute detailed alignment information from the given information.

    Alignment details are based off the differences between the alignment
    strings.
    This method returns a tuple containing:
        - The Start Position in the reference of the alignment.
        - The Template Length of this alignment.
        - The Cigar representing this alignment.
        - The Phred-Scaled quality score of this alignment.
    Where:
        - The Start Position is the 0-based, inclusive position in the reference
          at which this alignment begins.
        - The End Position is the 0-based, inclusive position in the reference
          at which this alignment ends.
        - The Template Length is the number of bases accounted by this alignment
          with respect to the reference.
        - The Cigar is a list of tuples: (CIGAR_ELEMENT, COUNT) where each
          CIGAR_ELEMENT is defined in pysam.
        - The Phred-Scaled quality score is defined by the following formula:
            -10 log_10((# mismatches + # insertions + # deletions)/target_length)
            However, because of how Tesserae works, we ignore leading and trailing deletions.
    """

    no_alignment = DetailedAlignmentInfo(-1, -1, 0, (), MIN_ALIGNMENT_PL)

    # Do some basic checks here:
    if len(query_alignment_string) == 0:
        return no_alignment

    start_index = get_start_index_from_alignment_start_string(target_alignment_string)

    # Now that we know where in the reference this target begins, we can start
    # to loop through both alignment strings at the same time.
    # In this loop we will:
    #     construct a cigar string
    #     determine counts for alignment quality score
    #     determine template length

    num_errors = 0

    cigar = []
    current_cigar_element = None
    current_cigar_element_count = 0

    num_leading_deletions = 0
    num_trailing_deletions = 0

    num_query_bases_used = 0

    have_only_seen_deletions = True

    for query_base, target_base in zip(
        query_alignment_string[start_index:].upper(), target_alignment_string[start_index:].upper()
    ):

        # The Tesserae2 / Mosaic Alignment algorithm can only produce "-" or
        # <BASE> for any position (other than blanks / spaces).  Therefore we
        # only have to check the following 4 cases:
        if query_base == "-":
            # We have an insertion relative to the reference:
            num_errors += 1
            cigar_element = pysam.CINS
            have_only_seen_deletions = False
            num_trailing_deletions = 0

        elif query_base == target_base:
            # Bases match:
            # We use CMATCH here because that cigar operator accounts for
            # BOTH matches and mismatches.
            cigar_element = pysam.CEQUAL
            have_only_seen_deletions = False
            num_trailing_deletions = 0
            num_query_bases_used += 1

        elif target_base == "-":
            # We have a deletion relative to the reference:
            num_errors += 1
            cigar_element = pysam.CDEL
            if have_only_seen_deletions:
                num_leading_deletions += 1
            num_trailing_deletions += 1
            num_query_bases_used += 1

        else:
            # We have a mismatch relative to the reference:
            num_errors += 1
            # We use CMATCH here because that cigar operator accounts for
            # BOTH matches and mismatches.
            cigar_element = pysam.CMATCH
            have_only_seen_deletions = False
            num_trailing_deletions = 0
            num_query_bases_used += 1

        # Accumulate our cigar elements:
        if cigar_element != current_cigar_element:
            if current_cigar_element is not None:
                cigar.append((current_cigar_element, current_cigar_element_count))

            current_cigar_element = cigar_element
            current_cigar_element_count = 1

        else:
            current_cigar_element_count += 1

    # Add the last remaining cigar element to our list:
    cigar.append((current_cigar_element, current_cigar_element_count))

    # Adjust cigar for leading / trailing deletions:
    if num_leading_deletions != 0:
        cigar = cigar[1:]
    if num_trailing_deletions != 0:
        cigar = cigar[:-1]

    # Our template length is the number of bases accounted by this alignment
    # with respect to the reference.
    # We add one because the end index is inclusive.
    template_length = target_end_index - target_start_index + 1

    # Compute end index (subtract 1 because of inclusive coordinates):
    end_index = start_index + num_query_bases_used - num_leading_deletions - num_trailing_deletions - 1

    # Make sure we have something to work with:
    if end_index < start_index:
        return no_alignment

    num_errors -= num_leading_deletions
    num_errors -= num_trailing_deletions

    # Compute PL score:
    aligned_qual_pl = get_qual_pl(num_errors, template_length)

    # Return our detailed info.
    return DetailedAlignmentInfo(start_index, end_index, template_length, tuple(cigar), aligned_qual_pl)


def create_pretty_alignment_string(query_alignment_info, target_aligned_info):
    """
    Create a pretty string containing all the information comparing the target alignments to the query alignment.

    :param query_alignment_info: A TesseraeAlignmentResult object representing the query string.
    :param target_aligned_info: A list[TesseraeAlignmentResult] representing the alignments of all targets to the
    query.
    :return: A string representing the alignment of all given targets against the given query.
    """

    out_string_list = ["Alignment:"]

    pos_label = "Position"

    # Get the width of the name field:
    name_width = max(len(query_alignment_info.seq_name), max([len(a.seq_name) for a in target_aligned_info]))

    # Get width of the numbering of the targets:
    tgt_num_width = math.ceil(math.log10(len(target_aligned_info)))

    # To deal with identifying the length of the string, we add numbers above the alignment representing the position
    # of each base and add start with that string.
    align_len = len(query_alignment_info.alignment_string)
    num_places = math.floor(math.log10(align_len))
    have_printed_pos_label = False
    for exponent in range(num_places, -1, -1):
        s = np.full(align_len, " ", dtype=str)
        ten_pow = int(math.pow(10, exponent))
        for i in range(0, align_len, min(ten_pow, 10)):
            s[i] = (i / ten_pow) % 10
        s[-1] = ((align_len - 1) / ten_pow) % 10

        if not have_printed_pos_label:
            out_string_list.append(f"{pos_label:>{tgt_num_width + 2 + name_width}} : " + "".join(s))
            have_printed_pos_label = True
        else:
            out_string_list.append(((tgt_num_width + 2 + name_width+3) * " ") + "".join(s))

    # Now that we have our numbered strings, we can add in the query sequence and the alignment string:
    out_string_list.append(f"{(tgt_num_width + 2) * ' '}{query_alignment_info.seq_name:>{name_width}} : "
                           f"{query_alignment_info.alignment_string}")

    # Default the alignments to bases being equal (we're optimists):
    s = np.full(align_len, " ", dtype=str)

    for t in target_aligned_info:
        # The starting position is the index of the first character in the alignment string of the first target that
        # is not a space:
        pos = get_start_index_from_alignment_start_string(t.alignment_string)

        # We take advantage of the fact that there will be only ONE alignment to each position and they are ordered:
        for t_base in t.alignment_string.lower()[pos:]:
            q_base = query_alignment_info.alignment_string[pos].lower()
            if q_base == "-":
                s[pos] = "~"
            elif t_base == "-":
                s[pos] = "^"
            elif t_base == q_base:
                s[pos] = "|"
            pos += 1
    out_string_list.append(((tgt_num_width + 2 + name_width+3) * " ") + "".join(s))

    # Now append all the targets:
    for i, target in enumerate(target_aligned_info):
        padding = " " * (align_len - len(target.alignment_string))
        out_string_list.append(f"{i:>{tgt_num_width}}: {target.seq_name:>{name_width}} : "
                               f"{target.alignment_string}{padding}")

    return "\n".join(out_string_list)


def get_start_index_from_alignment_start_string(target_alignment_string):
    """Get the alignment start index from the given alignment string.

    We know that the alignment strings will always start with spaces until
    the region that actually aligns to the reference, so we count the number
    of spaces in the target_alignment_string to get our starting alignment
    position.
    """

    try:
        start_index = next(i for i, v in enumerate(target_alignment_string) if v != " ")
    except StopIteration:
        # If we have a string of all spaces, we just need to skip it:
        start_index = len(target_alignment_string) - 1

    return start_index


def process_raw_results(raw_results_list, minqual, minbases):
    """
    Read through the raw results, computing additional statistics and filtering reads we shouldn't consider.

    :param raw_results_list: list of TesseraeAlignmentResults.
    :param minqual: minimum quality alignment to include (Phred scaled likelihood)
    :param minbases: minimum number of bases for an alignment to be retained.
    :return: A list of ProcessedAlignmentResult
    """

    # The read should be the first entry in the alignment:
    read_result = raw_results_list[0]  # TesseraeAlignmentResult

    last_target_name = None
    last_seq = None
    last_start_indx = None
    last_end_index = None

    added_last = False

    results = []
    for target_name, sequence, target_start_index, target_end_index in raw_results_list[1:]:

        # Merge adjacent sequences together!
        # Rules for merging:
        # Same sequence identifier
        # last sequence ends at base just before current sequence starts
        # last sequence does not end in a deletion ("-")
        # current sequence does not start with a deletion ("-")
        if target_name == last_target_name \
            and target_start_index == (last_end_index + 1) \
                and last_seq[-1] != "-" and sequence[0] != "-":

            LOGGER.info("Merging adjacent alignments: %s: [%d - %d] + [%d - %d]", target_name,
                        last_start_indx, last_end_index, target_start_index, target_end_index)

            # Now we can merge the sequences together.
            # We append the sequences together trimming whitespace from appropriate parts,
            # then we overwrite the start offset:
            sequence = last_seq + sequence.lstrip()
            target_start_index = last_start_indx

            LOGGER.debug("    Merged sequence: %s", sequence)
            LOGGER.debug("    Merged target start index: %d", target_start_index)

            # We must also remove the most recent entry in our results list, since it corresponds to the
            # entry for the "first half" of this alignment:
            # TODO: Fix this loop to not duplicate the alignment work.
            if added_last and results[-1].seq_name == target_name:
                results.pop()

        LOGGER.debug("Seq: %s", target_name)
        detailed_alignment_info = compute_detailed_alignment_info(
            read_result.alignment_string, sequence, target_start_index, target_end_index
        )

        if detailed_alignment_info.qual_pl > minqual and detailed_alignment_info.template_length >= minbases:

            processed_result = ProcessedAlignmentResult(
                target_name, sequence, target_start_index, target_end_index,
                detailed_alignment_info.start_pos, detailed_alignment_info.end_pos,
                detailed_alignment_info.template_length, detailed_alignment_info.cigar,
                detailed_alignment_info.qual_pl)

            LOGGER.debug("Adding in sequence: %s(q=%d-%d, tgt=%d-%d)", processed_result.seq_name,
                         processed_result.read_start_pos, processed_result.read_end_pos,
                         processed_result.target_start_index, processed_result.target_end_index)
            results.append(processed_result)
            added_last = True
        else:
            if LOGGER.isEnabledFor(logging.DEBUG):
                reason_string = ""
                if detailed_alignment_info.template_length < minbases:
                    reason_string = "[min aligned base: " + \
                                    str(detailed_alignment_info.template_length) + \
                                    " < " + \
                                    str(minbases) + "]"
                if detailed_alignment_info.qual_pl <= minqual:
                    sep = ""
                    if len(reason_string) > 0:
                        sep = " ,"
                    reason_string = reason_string + sep + "[quality - " + str(detailed_alignment_info.qual_pl) + \
                                    " <= " + str(minqual) + "]"
                LOGGER.debug("Target does not pass threshold: %s: %s (%s)",
                             target_name, reason_string, detailed_alignment_info)
            added_last = False

        # Track our info to the next iteration:
        last_target_name = target_name
        last_seq = sequence
        last_start_indx = target_start_index
        last_end_index = target_end_index

    # Since the alignments come out in the order they appear in the read sequence,
    # we should NOT sort them.
    # This will make sequence consolidation and order detection easier.
    return results


def remove_overlapping_bounded_seqs(bounded_seq_tuple_list):
    """
    Filters the given bounded sequence tuple list such that no sequences overlap.

    If two bounded regions overlap, one is selected as follows:
        start seg alignment quality (high to low)
        end seg alignment quality (high to low)
        start seg name (ABC Order)
        end seg name (ABC Order)
        Order in the given tuple list.

    :param bounded_seq_tuple_list: A list of tuples, with each tuple containing exactly two
    ProcessedAlignmentResult and representing an alignment to be excised from the parent read.
    :return: A set of tuples, with each tuple containing exactly two
    ProcessedAlignmentResult and representing an alignment to be excised from the parent read.
    Where no tuple overlaps the bounds of another tuple.
    """

    filtered_regions = []
    filtered_region_set = set()

    overlap_list_list = []
    overlapping_region_map = set()

    # 1 - find all overlapping tuples.
    for i, (result_start, result_end) in enumerate(bounded_seq_tuple_list):

        overlaps = []

        # Adjust for forward / reverse direction:
        if result_start.read_start_pos < result_start.read_end_pos:
            start = result_start.read_start_pos
            end = result_end.read_end_pos
        else:
            start = result_end.read_end_pos
            end = result_start.read_start_pos

        for other_result_start, other_result_end in bounded_seq_tuple_list[i+1:]:

            # Adjust for forward / reverse direction:
            if other_result_start.read_start_pos < other_result_end.read_end_pos:
                other_start = other_result_start.read_start_pos
                other_end = other_result_end.read_end_pos
            else:
                other_start = other_result_end.read_end_pos
                other_end = other_result_start.read_start_pos

            # Other fragment overlaps the start
            # Other fragment overlaps the end
            # current fragment overlaps the other's start
            # current fragment overlaps the other's end
            # current fragment is the same as the other's end
            if ((other_result_start, other_result_end) not in overlapping_region_map) and \
               (other_start < start < other_end) or (other_start < end < other_end) or \
               (start < other_start < end) or (start < other_end < end) or \
               (start == other_start and end == other_end):
                overlaps.append((other_result_start, other_result_end))

        # Do our best not to account for regions more than once:
        if (result_start, result_end) in overlapping_region_map:
            continue
        elif len(overlaps) == 0:
            filtered_regions.append((result_start, result_end))
            filtered_region_set.add((result_start, result_end))
        else:
            overlaps.insert(0, (result_start, result_end))
            overlap_list_list.append(overlaps)

            for r in overlaps:
                overlapping_region_map.add(r)

    LOGGER.debug("Overlap list lengths: %s", [len(x) for x in overlap_list_list])
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Overlap lists:")
        for i, tuple_list in enumerate(overlap_list_list):
            LOGGER.debug("List %d:", i)
            for t in tuple_list:
                LOGGER.debug("    (%s -> %s)", t[0].seq_name, t[1].seq_name)

    # Now we have a list of overlapping regions and non-overlapping regions.
    # We must resolve the overlaps, and chose only one alignment from each of the overlap lists:
    for overlap_list in overlap_list_list:
        current_best_region = overlap_list[0]
        best_read_region_len = current_best_region[1].read_end_pos + current_best_region[0].read_start_pos
        best_avg_qual = (current_best_region[0].overall_quality + current_best_region[1].overall_quality) / 2
        best_template_len = current_best_region[0].template_length + current_best_region[1].template_length

        for overlap_region in overlap_list[1:]:
            read_region_len = current_best_region[1].read_end_pos + current_best_region[0].read_start_pos
            avg_qual = (current_best_region[0].overall_quality + current_best_region[1].overall_quality) / 2
            template_len = current_best_region[0].template_length + current_best_region[1].template_length

            new_best = False
            if read_region_len > best_read_region_len:
                new_best = True
            elif read_region_len == best_read_region_len:
                if avg_qual > best_avg_qual:
                    new_best = True
                elif avg_qual == best_avg_qual:
                    if template_len > best_template_len:
                        new_best = True

            if new_best:
                best_read_region_len = read_region_len
                best_avg_qual = avg_qual
                best_template_len = template_len
                current_best_region = overlap_region

        # We have our best.  Now we add it:
        filtered_regions.append(current_best_region)
        LOGGER.debug(" Best region: %s", current_best_region)

    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Best region list:")
        for r in filtered_regions:
            LOGGER.debug("    %s", r)

    return filtered_regions


def find_all_bounded_seqs_in_alignment_results(processed_results, seq_boundaries):
    """
    Search the given processed results for sequences with the given seq_boundaries.

    Accounts for reverse-complemented read sections as well as the given (assumed forward) direction
    (as specified by the RC_READ_NAME_IDENTIFIER suffix).

    :param processed_results: A list of ProcessedAlignmentResult containing alignments to search.
    :param seq_boundaries: List of tuples, with each tuple containing exactly two sequence names.
    :return: A list of tuples, with each tuple containing exactly two ProcessedAlignmentResult and representing
    an alignment to be excised from the parent read.
    """

    hits = []

    start_seqs = [b[0] for b in seq_boundaries]
    end_seqs = [b[1] for b in seq_boundaries]

    # Go through the results and get the position of the first start sequence.
    # Then get the first position of corresponding end sequence.
    # Continue until you get to the end of the list.
    is_open_segment = False
    cur_start_result = None
    cur_indx = 0
    for i, result in enumerate(processed_results):
        if is_open_segment:
            if (result.seq_name in end_seqs) and ((cur_start_result.seq_name, result.seq_name) in seq_boundaries):
                hits.append((cur_start_result, result))
                is_open_segment = False
                LOGGER.debug("Found bounded subread: [%d - %d]", cur_indx, i)
                LOGGER.debug("    %s", cur_start_result)
                LOGGER.debug("    %s", result)

        elif result.seq_name in start_seqs:
            cur_start_result = result
            cur_indx = i
            is_open_segment = True

    num_forward_hits = len(hits)
    LOGGER.debug("Found %d forward direction segments for extraction.", num_forward_hits)

    # Now we account for the RC versions:
    is_open_segment = False
    cur_start_result = None
    cur_start_name_f = None
    cur_indx = 0
    for i, result in enumerate(processed_results[::-1]):

        # Correct the index for our traversal direction:
        i = len(processed_results) - i - 1

        # Get equivalent Forward sequence name:
        seg_name_f = result.seq_name[:-len(RC_READ_NAME_IDENTIFIER)]

        if is_open_segment:
            if (seg_name_f in end_seqs) and \
                    ((cur_start_name_f, seg_name_f) in seq_boundaries):
                hits.append((cur_start_result, result))
                is_open_segment = False
                LOGGER.debug("Found bounded subread: [%d - %d]", cur_indx, i)
                LOGGER.debug("    %s", cur_start_result)
                LOGGER.debug("    %s", result)

        elif seg_name_f in start_seqs:
            cur_start_result = result
            cur_start_name_f = seg_name_f
            cur_indx = i
            is_open_segment = True

    LOGGER.debug("Found %d reverse complement direction segments for extraction.", len(hits) - num_forward_hits)

    return hits


def dump_seq_map(seq_map, name='Reads'):
    """Dumps the given sequence map to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug('%s:', name)
        for e in seq_map.items():
            LOGGER.debug('    %s -> %s', e[0], e[1])


def dump_sequence_boundaries(seq_boundaries):
    """Dumps the given sequence list to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug('Sequence Boundaries:')
        for i, (s, e) in enumerate(seq_boundaries):
            LOGGER.debug('    % 2d: %s : %s', i+1, s, e)


def dump_results(results_tuple_list):
    """Dumps the results tuple to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug('Results:')
        for r in results_tuple_list:
            LOGGER.debug('    %s', str(r))


def assert_valid_sequence_boundaries(seq_boundaries, known_segment_names_to_seq_dict):
    """
    Validate that each boundary list in the given seq_boundaries contains only known sequences
    in the given known_segment_names_to_seq_dict.

    :param seq_boundaries: List of tuples, with each tuple containing exactly two sequence names.
    :param known_segment_names_to_seq_dict: dictionary of sequence name -> sequence upon which to validate seq_list.
    :return: None.  Raises a ValueError if the given seq_list is not valid.
    """

    for i, (s, e) in enumerate(seq_boundaries):
        if s not in known_segment_names_to_seq_dict:
            raise ValueError(
                f"Boundary {i+1} start is not in the known sequence list: {s}.  "
                "Sequence names must correspond to a sequence in the given Segments file."
            )
        elif e not in known_segment_names_to_seq_dict:
            raise ValueError(
                f"Boundary {i+1} end is not in the known sequence list: {e}.  "
                "Sequence names must correspond to a sequence in the given Segments file."
            )


def create_alignment_targets(known_segments_to_seq_dict, sequence_boundaries):
    """
    Create Tesserae Sequence objects to align known segments with Tesserae.

    Adds both the given direction and the reverse-complemented direction so that the sequences can be aligned in either
    direction.  This is significantly faster than running the whole alignment twice.

    :param known_segments_to_seq_dict: dict of sequence_name -> sequence_bases
    :param sequence_boundaries: List of tuples, with each tuple containing exactly two sequence names.
    :return: list(tesserae.Sequence) of all known sequences that are not UNKNOWN / BARCODE.
    """

    alignment_target_seqs = []

    for bound_seqs in sequence_boundaries:
        for s in bound_seqs:
            alignment_target_seqs.append(Sequence(s, known_segments_to_seq_dict[s]))
            alignment_target_seqs.append(
                Sequence(s + RC_READ_NAME_IDENTIFIER, reverse_complement(known_segments_to_seq_dict[s]))
            )

    return alignment_target_seqs


def ingest_fastx_file(file_path):
    """Ingest the contents of a FASTA/FASTQ file and return two dictionaries
    onc
        1: Mapping from read name to read sequence
        2: Mapping from template name to read name

    The `template name` is simply the word 'template' with the
    order in which a given read occurs in the file
    (e.g. 'template0' or 'template10').
    """

    t_num = 0
    _read_to_sequence_dict = OrderedDict()
    _template_to_read_name_dict = dict()
    with pysam.FastxFile(file_path) as file_handle:
        for entry in file_handle:
            _read_to_sequence_dict[entry.name] = entry.sequence
            _template_to_read_name_dict[f"template{t_num}"] = entry.name

            t_num += 1

    return _read_to_sequence_dict, _template_to_read_name_dict


def ingest_sequence_boundaries(file_path):
    """Ingest the contents of a boundaries file.

    A boundaries file is a plain text file with two comma separated sequence names per line.
    These sequence names should correspond to the names of sequences in a FASTX
    file that is given as an argument to Cartographer.

    Sequence names will have all whitespace at the start and end stripped off (after splitting by ",")

    The validity of sequence names is not checked here."""

    seq_boundaries = []
    with open(file_path, 'r') as f:
        for i, line in enumerate(f.readlines()):
            bounds = line.strip().split(",")
            if len(bounds) != 2:
                raise ValueError(
                    "Each line in the boundaries file must contain exactly two sequences.  "
                    f"On line {i} found {len(bounds)}",
                )
            seq_boundaries.append((bounds[0].strip(), bounds[1].strip()))

    return seq_boundaries


################################################################################


def extract_read_sections(args):
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

    alignment_algorithm = AlignmentAlgorithm.TESSERAE
    if args.aligner:
        try:
            alignment_algorithm = AlignmentAlgorithm[args.aligner]
        except KeyError:
            LOGGER.error("Error: You must provide a valid alignment algorithm.  Options are: %s",
                         ", ".join([e.name for e in AlignmentAlgorithm]))
            sys.exit(1)
    LOGGER.info("Alignment Algorithm: %s", alignment_algorithm.name)

    LOGGER.info("Ingesting boundaries from %s ...", args.boundaries)
    sequence_boundaries = ingest_sequence_boundaries(args.boundaries)
    LOGGER.info("Ingested %d sequence boundaries.", len(sequence_boundaries))
    dump_sequence_boundaries(sequence_boundaries)

    LOGGER.info("Ingesting known segments from %s ...", args.segments)
    known_segment_names_to_seq_dict, _ = ingest_fastx_file(args.segments)
    LOGGER.info("Ingested %d known segments.", len(known_segment_names_to_seq_dict))
    dump_seq_map(known_segment_names_to_seq_dict, "Known Segments")

    LOGGER.info("Validating given sequence boundaries...")
    assert_valid_sequence_boundaries(sequence_boundaries, known_segment_names_to_seq_dict)
    LOGGER.debug("Boundary list is valid.")

    # Create an ordered dict of our ordered sequences:
    alignment_target_seqs = create_alignment_targets(known_segment_names_to_seq_dict, sequence_boundaries)

    # A couple of spacing variables for nice looking logs:
    spacing_one = " " * 4

    # Open all our files here so they'll be automatically closed:
    with open(args.outfile, 'w') as out_file, \
            open(args.rejected_outfile, 'w') as rejected_out_file, \
            open(args.raw_marker_alignments, 'w') as raw_marker_alignments_file, \
            open(args.initial_section_alignments, 'w') as initial_section_alignments, \
            open(args.final_section_alignments, 'w') as final_section_alignments_file:

        num_sequences_extracted = 0
        num_reads_with_sub_sequences = 0
        num_forward_subsequences_extracted = 0
        num_rc_subsequences_extracted = 0
        num_rejected = 0
        LOGGER.info("Processing reads...")

        with ReadFile(args.reads) as reads_file:
            num_reads = 0
            for read_num, read_sequence in enumerate(reads_file.get_reads()):
                num_reads += 1

                LOGGER.debug(
                    "%sProcessing read %d: %s (len=%d)",
                    spacing_one,
                    read_num,
                    read_sequence.name,
                    len(read_sequence.sequence)
                )

                if args.max_read_length and len(read_sequence.sequence) > args.max_read_length:
                    LOGGER.warning("Ignoring read %d - %s: Length too long: %d > %d", read_num, read_sequence.name,
                                   len(read_sequence.sequence), args.max_read_length)
                    rejected_out_file.write(f">{read_sequence.name}\n")
                    rejected_out_file.write(f">{read_sequence.sequence}\n")
                    num_rejected += 1
                    continue

                bounded_seq_tuple_list = []

                LOGGER.info("%sInitial alignment of known segments ...", spacing_one)
                segment_alignment_results, query_result = align_sequences(
                    read_sequence, alignment_target_seqs, args.minqual, args.minbases, p_rec=P_REC_KNOWN,
                    alignment_type=alignment_algorithm, threads=num_threads
                )

                if len(segment_alignment_results) != 0:
                    LOGGER.debug("%sDetermining if expected sequence appears in processed alignment results...",
                                 spacing_one)

                    # Create a position map so we can get real read positions from teh alignment string positions:
                    alignment_read_pos_map = IdentityMap()
                    if alignment_algorithm != AlignmentAlgorithm.BWA_MEM and \
                            alignment_algorithm != AlignmentAlgorithm.BWA_ALN:
                        # We only need to adjust the positions we're not using the BWA aligners:
                        alignment_read_pos_map = create_alignment_to_base_map(query_result.alignment_string)

                    # Write our raw marker alignments:
                    write_marker_alignments_to_output_file(
                        raw_marker_alignments_file,
                        read_sequence.name,
                        segment_alignment_results,
                        alignment_read_pos_map
                    )

                    bounded_seq_tuple_list = find_all_bounded_seqs_in_alignment_results(
                        segment_alignment_results, sequence_boundaries
                    )

                    # Write our initial section alignments:
                    write_section_alignments_to_output_file(
                        initial_section_alignments,
                        read_sequence.name,
                        bounded_seq_tuple_list,
                        alignment_read_pos_map
                    )

                LOGGER.info("%sExpected ordered sequences occur %d time(s).", spacing_one, len(bounded_seq_tuple_list))

                if len(bounded_seq_tuple_list) == 0:
                    LOGGER.info("%sThis read has no complete matches for the given known ordered sequences.",
                                spacing_one)
                else:
                    # We now have a list of ALL matching sequences.
                    # We need to filter it so that none of the matches overlap
                    # (that would double-count and would also be silly).
                    filtered_bounded_seq_tuple_list = remove_overlapping_bounded_seqs(bounded_seq_tuple_list)

                    # Write our initial section alignments:
                    write_section_alignments_to_output_file(
                        final_section_alignments_file,
                        read_sequence.name,
                        bounded_seq_tuple_list,
                        alignment_read_pos_map
                    )

                    # Write out our new subsequences to the output fasta file:
                    write_sub_sequences(read_sequence, filtered_bounded_seq_tuple_list, alignment_read_pos_map, out_file)

                    # Track subsequence statistics:
                    num_forward_subsequences_extracted += sum(
                        not b[0].seq_name.endswith(RC_READ_NAME_IDENTIFIER) for b in filtered_bounded_seq_tuple_list
                    )
                    num_rc_subsequences_extracted += sum(
                        b[0].seq_name.endswith(RC_READ_NAME_IDENTIFIER) for b in filtered_bounded_seq_tuple_list
                    )
                    num_sequences_extracted_this_read = len(filtered_bounded_seq_tuple_list)
                    num_sequences_extracted += num_sequences_extracted_this_read
                    num_reads_with_sub_sequences += 1

            LOGGER.info("Processed %d reads.", num_reads)
            if args.max_read_length:
                LOGGER.info("Rejected %d reads.", num_rejected)
            LOGGER.info("# Reads containing sub-sequences: %d", num_reads_with_sub_sequences)
            LOGGER.info("# forward direction sub-sequences extracted: %d", num_forward_subsequences_extracted)
            LOGGER.info("# reverse-complemented direction sub-sequences extracted: %d", num_rc_subsequences_extracted)
            LOGGER.info("Total # sub-sequences extracted: %d", num_sequences_extracted)


def write_marker_alignments_to_output_file(out_file, read_name, alignments, align_pos_read_pos_map):
    """
    Writes the given marker alignments to the given out_file.
    Writes one line: the read_name, followed by a basic string representation of each alignment separated by tabs.

    If len(alignments) == 0, then does not write anything.

    :param out_file: An open File object to which to write the data.
    :param read_name: The name of the read to which the given alignments belong.
    :param alignments: list(ProcessedAlignmentResult) representing alignments to write to the file.
    :param align_pos_read_pos_map: Map from alignment string position to read position.
    """

    if len(alignments) > 0:
        out_file.write(read_name)
        for a in alignments:
            LOGGER.debug(f"a.seq_name = {a.seq_name}")
            LOGGER.debug(f"a.read_start_pos = {a.read_start_pos}")
            LOGGER.debug(f"a.read_end_pos = {a.read_end_pos}")
            LOGGER.debug(f"a.overall_quality = {a.overall_quality}")
            out_file.write(f"\t{a.seq_name}:"
                           f"{align_pos_read_pos_map[a.read_start_pos]}"
                           f"-{align_pos_read_pos_map[a.read_end_pos]}"
                           f"@{a.overall_quality}")
        out_file.write("\n")


def write_section_alignments_to_output_file(out_file, read_name, section_tuples, align_pos_read_pos_map):
    """
    Writes the given marker alignments to the given out_file.
    Writes one line: the read_name, followed by a basic string representation of each alignment separated by tabs.

    If len(section_tuples) == 0, then does not write anything.

    :param out_file: An open File object to which to write the data.
    :param read_name: The name of the read to which the given alignments belong.
    :param section_tuples: A list of tuples, with each tuple containing exactly two ProcessedAlignmentResult and
    representing an alignment to be excised from the parent read.
    :param align_pos_read_pos_map: Map from alignment string position to read position.
    """
    if len(section_tuples) != 0:
        out_file.write(read_name)
        for s1, s2 in section_tuples:
            out_file.write("\t")
            out_file.write(f"[{s1.seq_name}:"
                           f"{align_pos_read_pos_map[s1.read_start_pos]}"
                           f"-{align_pos_read_pos_map[s1.read_end_pos]}"
                           f"@{s1.overall_quality}")
            out_file.write("<>")
            out_file.write(f"{s2.seq_name}:"
                           f"{align_pos_read_pos_map[s2.read_start_pos]}"
                           f"-{align_pos_read_pos_map[s2.read_end_pos]}"
                           f"@{s2.overall_quality}]")
        out_file.write("\n")


def create_alignment_to_base_map(alignment_string):
    """
    Create a dictionary mapping positions in the given alignment string back to the original string it came from.

    Assumes that every base in the original string appears in the alignment string and the only additional characters
    are '-' to represent insertions.

    :param alignment_string: Alignment string containing original bases or '-' (to indicate insertions).
    :return: A dict() mapping positions in the given alignment string back to the original string.
    """

    position_map = {i: i for i in range(len(alignment_string))}

    o_pos = 0
    for i, c in enumerate(alignment_string):
        position_map[i] = o_pos
        if not c == "-":
            o_pos += 1

    return position_map


def write_sub_sequences(read_sequence, bounded_seq_tuple_list, alignment_read_pos_map, out_fasta_file):
    """
    Write out the sections of the given read_sequence as bounded by bounded_seq_list to the given out_file in
    FASTA format.
    :param read_sequence: A tesserae.Sequence object representing the parent read to the given bounded_seq_list.
    :param bounded_seq_tuple_list: A list of tuples, with each tuple containing exactly two ProcessedAlignmentResult
    and representing an alignment to be excised from the parent read.
    :param alignment_read_pos_map: Dictionary mapping positions in the given alignment string back to the original
    string.  This is required to decode all positions from alignment string positions to original read positions.
    :param out_fasta_file: An open file object to which to write results.
    :return: None
    """

    for start_alignment, end_alignment in bounded_seq_tuple_list:

        # Do some math here to account for how we create the tuple list:
        # And adjust for reverse complements:
        start_coord = alignment_read_pos_map[start_alignment.read_start_pos]
        end_coord = alignment_read_pos_map[end_alignment.read_end_pos]
        start_name = start_alignment.seq_name
        end_name = end_alignment.seq_name
        if start_alignment.seq_name.endswith(RC_READ_NAME_IDENTIFIER):
            start_coord = end_alignment.read_end_pos
            end_coord = start_alignment.read_start_pos
            # Clean up the disambiguators (md5sums) from the RC names:
            decorator_length = len(RC_READ_NAME_IDENTIFIER) - len(RC_READ_NAME_IDENTIFIER_BASE)
            start_name = start_alignment.seq_name[:-decorator_length]
            end_name = end_alignment.seq_name[:-decorator_length]

        out_fasta_file.write(f">{read_sequence.name}_{start_coord}-{end_coord}_"
                             f"{start_name}-{end_name}\n")

        # Quick bounds check on the coords for the read sequence:
        if start_coord < end_coord:
            out_fasta_file.write(f"{read_sequence.sequence[start_coord:end_coord+1]}\n")
        else:
            out_fasta_file.write(f"{read_sequence.sequence[end_coord:start_coord+1]}\n")


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


def create_alignment_with_bwa_mem(read_sequence, target_sequences, minqual, minbases, threads=1, log_spacing="    "):
    """
    Perform a BWA MEM 2 alignment on

    :param read_sequence: tesserae Sequence object against which to align all sequences in target_sequences
    :param target_sequences: ordered list of tessereae Sequence objects to align.
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :param threads: number of threads to use when aligning.
    :param log_spacing: Spacing to precede any log statements.
    :return: A list of ProcessedAlignmentResult objects.
    """

    out_file_name = "tmp.sam"

    # Write sequences to tmp fasta file:
    _, ref_file = tempfile.mkstemp()
    _, seq_file = tempfile.mkstemp()
    try:
        LOGGER.debug("%sCreating tmp \"reference\" fasta file: %s", log_spacing, ref_file)
        with open(ref_file, "w", encoding="ascii") as tmp:
            tmp.write(f">{read_sequence.name}\n")
            tmp.write(f"{read_sequence.sequence}\n")

        bwa_index_args = ["/bwa-mem2-2.0pre2_x64-linux/bwa-mem2", "index", ref_file]
        LOGGER.debug(
            "%sRunning BWA Index on tmp file (%s): %s", log_spacing, ref_file, " ".join(bwa_index_args)
        )
        _ = subprocess.run(bwa_index_args, capture_output=True, check=True)

        LOGGER.debug("%sCreating tmp known segment \"read\" fasta file: %s", log_spacing, seq_file)
        seen_targets = set()
        with open(seq_file, "w", encoding="ascii") as tmp:
            for t in target_sequences:
                # We don't need to explicitly align reverse-complemented sequences here:
                if t.name.endswith(RC_READ_NAME_IDENTIFIER):
                    continue
                elif t.name not in seen_targets:
                    tmp.write(f">{t.name}\n")
                    tmp.write(f"{t.sequence}\n")
                    seen_targets.add(t.name)

        LOGGER.debug("Contents of tmp \"reference\" fasta file:")
        with open(ref_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        LOGGER.debug("Contents of known segment \"read\" fasta file:")
        with open(seq_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        bwa_mem_args = ["/bwa-mem2-2.0pre2_x64-linux/bwa-mem2", "mem",
                        "-a",                # Output all found alignments for single-end or unpaired paired-end reads.
                                             # These alignments will be flagged as secondary alignments.
                        "-S",                # skip mate rescue
                        "-P",                # skip pairing; mate rescue performed unless -S also in use
                        "-k8",               # minimum seed length
                        "-A", "1",           # Matching score.
                        "-B", "4",           # Mismatch penalty.
                                             #  The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}.
                        "-O", "6,6",         # Gap open penalty.
                        "-E", "1,1",         # gap extension penalty; a gap of size k cost '{-O} + {-E}*k'
                        "-L", "5,5",         # penalty for 5'- and 3'-end clipping
                        "-U", "17",          # penalty for an unpaired read pair
                        "-T", "30",          # minimum score to output
                        "-c", "1000",        # skip seeds with more than INT occurrences
                        "-t", str(threads),
                        "-o", out_file_name,
                        ref_file, seq_file]
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
        os.remove(seq_file)
        try:
            os.remove(out_file_name)
            pass
        except FileNotFoundError:
            # If the alignment failed, we won't necessarily have an output file:
            pass


def create_alignment_with_bwa_aln(read_sequence, target_sequences, minqual, minbases, threads=1, log_spacing="    "):
    """
    Perform a BWA ALN alignment on

    :param read_sequence: tesserae Sequence object against which to align all sequences in target_sequences
    :param target_sequences: ordered list of tessereae Sequence objects to align.
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :param threads: number of threads to use when aligning.
    :param log_spacing: Spacing to precede any log statements.
    :return: A list of TesseraeAlignmentResult objects.
    """

    out_sai_file = "tmp.sai"
    out_file_name = "tmp.sam"

    # Write sequences to tmp fasta file:
    _, ref_file = tempfile.mkstemp()
    _, seq_file = tempfile.mkstemp()
    try:
        LOGGER.debug("%sCreating tmp \"reference\" fasta file: %s", log_spacing, ref_file)
        with open(ref_file, "w", encoding="ascii") as tmp:
            tmp.write(f">{read_sequence.name}\n")
            tmp.write(f"{read_sequence.sequence}\n")

        bwa_index_args = ["/bwa/bwa", "index", ref_file]
        LOGGER.debug(
            "%sRunning BWA Index on tmp file (%s): %s", log_spacing, ref_file, " ".join(bwa_index_args)
        )
        _ = subprocess.run(bwa_index_args, capture_output=True, check=True)

        LOGGER.debug("%sCreating tmp known segment \"read\" fasta file: %s", log_spacing, seq_file)
        seen_targets = set()
        with open(seq_file, "w", encoding="ascii") as tmp:
            for t in target_sequences:
                # We don't need to explicitly align reverse-complemented sequences here:
                if t.name.endswith(RC_READ_NAME_IDENTIFIER):
                    continue
                elif t.name not in seen_targets:
                    tmp.write(f">{t.name}\n")
                    tmp.write(f"{t.sequence}\n")
                    seen_targets.add(t.name)

        LOGGER.debug("Contents of tmp \"reference\" fasta file:")
        with open(ref_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        LOGGER.debug("Contents of known segment \"read\" fasta file:")
        with open(seq_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        # Run the alignment:
        bwa_aln_args = [
            "/bwa/bwa", "aln",
            "-n 0.04",               # max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
            "-o 1",                  # maximum number or fraction of gap opens [1]
            "-e -1",                 # maximum number of gap extensions, -1 for disabling long gaps [-1]
            "-d 10",                 # maximum occurrences for extending a long deletion [10]
            "-k 2",                  # maximum differences in the seed [2]
            "-M 3",                  # mismatch penalty [3]
            "-O 11",                 # gap open penalty [11]
            "-E 4",                  # gap extension penalty [4]
            "-l 8",                  # seed length [32]
            "-i 1",                  # do not put an indel within INT bp towards the ends [5]
            "-R 30",                 # stop searching when there are >INT equally best hits
            "-q 0",                  # quality threshold for read trimming down to 35bp [0]
            "-t", str(threads),
            "-f", out_sai_file,
            ref_file,
            seq_file
        ]

        LOGGER.debug(
            "%sRunning BWA ALN on tmp file (%s): %s", log_spacing, ref_file, " ".join(bwa_aln_args)
        )
        subprocess.run(bwa_aln_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)

        # Convert the alignment to sam file:
        bwa_samse_args = [
            "/bwa/bwa", "samse",
            f"-n {len(target_sequences)}",  # Maximum number of alignments to output in the XA tag for
                                            # reads paired properly. If a read has more than INT hits, the
                                            # XA tag will not be written. [3]
            f"-f{out_file_name}",
            ref_file,
            out_sai_file,
            seq_file
        ]

        LOGGER.debug(
            f"{log_spacing}Running BWA SAMSE on tmp SAI file ({out_sai_file}): {' '.join(bwa_samse_args)}"
        )
        completed_process = subprocess.run(bwa_samse_args,
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)

        if LOGGER.isEnabledFor(logging.DEBUG):
            LOGGER.debug("BWA ALN output:")
            for l in completed_process.stdout.split("\n"):
                LOGGER.debug(l)

            with open(out_file_name, 'rb') as f:
                LOGGER.debug("=" * 80)
                LOGGER.debug("Raw BWA ALN Alignment:")
                for line in f.readlines():
                    LOGGER.debug("%s", line.decode("ascii").rstrip())
                LOGGER.debug("=" * 80)

        return get_processed_results_from_bwa_aln_file(out_file_name, minqual, minbases)

    except subprocess.CalledProcessError as e:
        LOGGER.error("Could not align with BWA ALN!")
        LOGGER.error("Stdout: %s", e.stdout.decode("utf-8"))
        LOGGER.error("Stderr: %s", e.stderr.decode("utf-8"))
        raise e

    finally:
        os.remove(ref_file)
        os.remove(seq_file)
        try:
            os.remove(out_sai_file)
            pass
        except FileNotFoundError:
            # If the alignment failed, we won't necessarily have an output file:
            pass
        try:
            os.remove(out_file_name)
            pass
        except FileNotFoundError:
            # If the alignment failed, we won't necessarily have an output file:
            pass


def create_raw_alignment_with_mosaic_aligner(read_sequence, target_sequences, log_spacing="    "):
    """
    Perform a MosaicAligner alignment on

    :param read_sequence: tesserae Sequence object against which to align all sequences in target_sequences
    :param target_sequences: ordered list of tessereae Sequence objects to align.
    :param log_spacing: Spacing to precede any log statements.
    :return: A list of TesseraeAlignmentResult objects.
    """

    out_file_name = "align.txt"

    # Write sequences to tmp fasta file:
    fd, seq_file = tempfile.mkstemp()
    LOGGER.debug("%sCreating tmp fasta file: %s", log_spacing, seq_file)
    seen_targets = set()
    try:
        with open(seq_file, "w", encoding="ascii") as tmp:
            tmp.write(f">{read_sequence.name}\n")
            tmp.write(f"{read_sequence.sequence}\n")

            for t in target_sequences:
                if t.name not in seen_targets:
                    tmp.write(f">{t.name}\n")
                    tmp.write(f"{t.sequence}\n")
                    seen_targets.add(t.name)

        LOGGER.debug("Contents of tmp fasta file:")
        with open(seq_file, "r") as f:
            for l in f.readlines():
                LOGGER.debug(l.rstrip())

        mosaic_aligner_args = ["/MosaicAligner/mosaic", "-breakOnLowerSeqPos", "-nt", "-seq", seq_file]
        LOGGER.debug(
            "%sRunning MosaicAligner on tmp file (%s): %s", log_spacing, seq_file, " ".join(mosaic_aligner_args)
        )
        _ = subprocess.run(mosaic_aligner_args, capture_output=False, check=True)

        if LOGGER.isEnabledFor(logging.DEBUG):
            with open(out_file_name, 'rb') as f:
                LOGGER.debug("=" * 80)
                LOGGER.debug("Raw MosaicAligner Alignment:")
                for line in f.readlines():
                    LOGGER.debug("%s", line.decode("ascii").rstrip())
                LOGGER.debug("=" * 80)

        return get_raw_results_from_mosaic_aligner_file(out_file_name)

    except subprocess.CalledProcessError as e:
        LOGGER.error("Could not align with MosaicAligner!")
        LOGGER.error("Stdout: %s", e.stdout.decode("utf-8"))
        LOGGER.error("Stderr: %s", e.stderr.decode("utf-8"))
        raise e

    finally:
        os.remove(seq_file)
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
            else:
                processed_results.append(p)

    # Sort by the order in which they appear in the read.
    # This is _VERY_IMPORTANT_ for finding the ordered regions in a following step.
    processed_results.sort(key=lambda x: x.read_start_pos)

    for r in processed_results:
        LOGGER.debug("    %s", r)

    return processed_results


def parse_cigar_string_to_tuples(cigar_string):
    """
    Parse the given cigar string into a list of cigar tuples.
    :param cigar_string: String containing CIGAR information.
    :return: Tuple of cigar tuples corresponding to the given CIGAR string.
    """
    # This is not very pythonic and I don't know how to fix it offhand...
    cigar_tuple_list = []
    buf = []
    for c in cigar_string:
        # Still in the number portion of the cigar element:
        v = ord(c)
        if 47 < v < 58:
            buf.append(c)
        else:
            # Now we're at the "string" portion:
            try:
                v = STRING_CIGAR_ELEMENT_MAP[c]
            except KeyError:
                raise KeyError(f"Error: No such cigar element: {c}")

            cigar_tuple_list.append((v, int("".join(buf))))
            buf = []

    return tuple(cigar_tuple_list)


def get_processed_results_from_bwa_aln_file(file_path, minqual, minbases):
    """
    Ingests the given file and creates raw alignments from it.
    This method is VERY similar to the BWA MEM equivalent, but with an epilog to handle XA tags for each read.

    :param file_path: Path to the sam file of a BWA ALN run.
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

            # Make a list of candidate read alignments so we can filter them all later:
            candidate_alignments = []

            # Note - must adjust ref start/end pos to align properly with conventions from other aligners
            #        (other aligner conventions: 1-based coordinates, inclusive end positions)
            p = ProcessedAlignmentResult(
                seq_name, bases, leading_clips, len(bases)-trailing_clips-1,
                int(read.reference_start), int(read.reference_end-1),
                template_length, tuple(read.cigartuples), qual_pl
            )
            candidate_alignments.append(p)

            # Get all XA alignments as well:
            candidate_alignments.extend(process_xa_tags(read, read_seqs))

            # Process all alignments with the same rules:
            for p in candidate_alignments:
                # Check against thresholds to make sure we should report the alignment:
                if p.overall_quality < minqual or p.template_length < minbases:
                    if p.overall_quality < minqual and p.template_length < minbases:
                        reason_string = f"qual too low ({p.overall_quality} < {minqual}) " \
                                        f"AND aligment too short ({p.template_length} < {minbases})"
                    elif p.template_length < minbases:
                        reason_string = f"aligment too short ({p.template_length} < {minbases})"
                    else:
                        reason_string = f"qual too low ({p.overall_quality} < {minqual})"

                    LOGGER.debug("Target does not pass threshold: %s: %s (%s)", p.seq_name, reason_string, p)
                else:
                    processed_results.append(p)

    # Sort by the order in which they appear in the read.
    # This is _VERY_IMPORTANT_ for finding the ordered regions in a following step.
    processed_results.sort(key=lambda x: x.read_start_pos)

    for r in processed_results:
        LOGGER.debug("    %s", r)

    return processed_results


def process_xa_tags(read, read_seqs):
    xa_processed_reads = []
    try:
        # split the XA tags by delimiters:
        alternate_alignments = read.get_tag("XA").split(';')
        for xa_alignment_string in alternate_alignments:
            if len(xa_alignment_string) == 0:
                continue

            LOGGER.debug("XA Alignment: %s", xa_alignment_string)

            # Parse the data into an aligned segment object:
            ref_name, pos, cigar, edit_dist = xa_alignment_string.split(',')
            pos = int(pos)
            edit_dist = int(edit_dist)
            seq_name = read.query_name

            xa_read = pysam.AlignedSegment()
            xa_read.query_name = seq_name
            xa_read.query_sequence = read_seqs[seq_name]
            xa_read.is_reverse = (pos < 0)
            xa_read.reference_start = pos if pos >= 0 else -pos
            xa_read.cigartuples = parse_cigar_string_to_tuples(cigar)
            xa_read.set_tag("NM", edit_dist, value_type='i')

            if pos < 0:
                seq_name = seq_name + RC_READ_NAME_IDENTIFIER

            template_length = xa_read.infer_read_length()
            qual_pl = get_qual_pl(edit_dist, template_length)

            # Get the leading and trailing clips so we can remove them from the aligned string:
            leading_clips = 0
            for e in xa_read.cigartuples:
                if e[0] == pysam.CSOFT_CLIP or e[0] == pysam.CHARD_CLIP:
                    leading_clips += e[1]
                else:
                    break

            trailing_clips = 0
            for e in xa_read.cigartuples[::-1]:
                if e[0] == pysam.CSOFT_CLIP or e[0] == pysam.CHARD_CLIP:
                    trailing_clips += e[1]
                else:
                    break

            # Note - must adjust ref start/end pos to align properly with conventions from other aligners
            #        (other aligner conventions: 1-based coordinates, inclusive end positions)
            p = ProcessedAlignmentResult(
                seq_name, xa_read.query_sequence, leading_clips, len(xa_read.query_sequence) - trailing_clips - 1,
                int(xa_read.reference_start), int(xa_read.reference_end-1),
                template_length, tuple(xa_read.cigartuples), qual_pl
            )
            LOGGER.debug(f"Adding: {p}")
            xa_processed_reads.append(p)

    except KeyError:
        return []

    return xa_processed_reads


def get_raw_results_from_mosaic_aligner_file(file_path="align.txt"):
    """
    Ingests the given file and creates raw alignments from it.
    :param file_path: Path to the results of a MosaicAligner run.
    :return: A list of TesseraeAlignmentResult objects.
    """

    raw_results = []

    with open(file_path, "r") as f:

        # The only result we need is the first one in the file:

        while not f.readline().startswith("Target"):
            pass

        # Now we're at our first alignment line (corresponding to the query line).
        # We can now iterate through the rest of the lines:
        line = f.readline().rstrip()
        while len(line) > 1:
            # Add in all our results:
            match = MOSAIC_ALIGNER_LINE_MATCH_PATTERN.match(line)
            if match:
                t = TesseraeAlignmentResult(
                        match.group(1),
                        match.group(4),
                        int(match.group(2)),
                        int(match.group(3)),
                    )
                raw_results.append(t)
            line = f.readline().rstrip()

    return raw_results


def align_sequences(read_sequence, target_sequences, minqual, minbases,
                    p_rec=DEFAULT_REC, alignment_type=AlignmentAlgorithm.TESSERAE,
                    threads=1, log_spacing="    "):
    """
    Perform a alignment on

    :param read_sequence: tesserae Sequence object against which to align all sequences in target_sequences
    :param target_sequences: ordered list of tessereae Sequence objects to align.
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :param p_rec: Prior probability of recombination event taking place.
    :param alignment_type: AlignmentAlgorithm object specifying which alignment algorithm to use.
    :param threads: Number of threads with which to run Tesserae alignment.
    :param log_spacing: Spacing to precede any log statements.
    :return: A tuple of results and the read sequence result (list(ProcessedAlignmentResult), TesseraeAlignmentResult)
    """

    # BWA MEM Alignments are sufficiently different that we have to handle them separately:
    if alignment_type == AlignmentAlgorithm.BWA_MEM:
        start_time = time.time()
        processed_results = create_alignment_with_bwa_mem(read_sequence, target_sequences, minqual, minbases,
                                                          threads=threads)
        end_time = time.time()
        LOGGER.info("%sCreated %d Alignments.  Alignment took %fs", log_spacing, len(processed_results),
                    end_time - start_time)

        query_result = TesseraeAlignmentResult(read_sequence.name, read_sequence.sequence, 0, len(read_sequence.sequence))
        return processed_results, query_result

    # BWA ALN Alignments are sufficiently different that we have to handle them separately:
    if alignment_type == AlignmentAlgorithm.BWA_ALN:
        start_time = time.time()
        processed_results = create_alignment_with_bwa_aln(read_sequence, target_sequences, minqual, minbases,
                                                          threads=threads)
        end_time = time.time()
        LOGGER.info("%sCreated %d Alignments.  Alignment took %fs", log_spacing, len(processed_results),
                    end_time - start_time)

        query_result = TesseraeAlignmentResult(read_sequence.name,
                                               read_sequence.sequence,
                                               0,
                                               len(read_sequence.sequence))
        return processed_results, query_result

    # Handle all other alignment types here:
    start_time = time.time()
    if alignment_type == AlignmentAlgorithm.MOSAIC_ALIGNER:
        raw_results = create_raw_alignment_with_mosaic_aligner(read_sequence, target_sequences)
    elif alignment_type == AlignmentAlgorithm.TESSERAE:
        raw_results = Tesserae(prho=p_rec, threads=threads).align(read_sequence, target_sequences)
    elif alignment_type == AlignmentAlgorithm.SMITH_WATERMAN:
        LOGGER.error("SMITH WATERMAN not implemented yet for structural alignment.")
        sys.exit(1)
    else:
        LOGGER.error("%s not implemented yet for structural alignment.", alignment_type.name)
        sys.exit(1)
    end_time = time.time()

    LOGGER.info("%sCreated %d Alignments.  Alignment took %fs", log_spacing, len(raw_results), end_time - start_time)
    dump_results(raw_results)

    query_result = raw_results[0]

    LOGGER.debug(create_pretty_alignment_string(query_result, raw_results[1:]))

    # Clean alignment results and calculate qualities:
    LOGGER.debug("%sProcessing raw results...", log_spacing)

    start_time = time.time()
    alignment_results = process_raw_results(raw_results, minqual, minbases)
    end_time = time.time()

    LOGGER.debug("%sProcessing results took %fs", log_spacing, end_time - start_time)
    LOGGER.debug("%sDumping processed results:", log_spacing)
    dump_results(alignment_results)

    if len(alignment_results) != 0:
        LOGGER.debug(create_pretty_alignment_string(query_result, alignment_results))

    return alignment_results, query_result


################################################################################


def main(raw_args):

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        description="Ingests three files: "
                    "one read file(FASTA/FASTQ/SAM/BAM) containing reads to be mapped, "
                    "one FASTA file containing known possible sequences that can occur in the read,"
                    "and a file containing read sequence boundaries.  This sequence boundaries file is"
                    "a plain text file with two comma separated sequence names per line.  The names should correspond"
                    "to the sequence names in the given sequence FASTA file.",
        usage="extract read sections bounded by given known sequences into a new fasta file",
    )

    align_required_args = parser.add_argument_group("required arguments")
    align_required_args.add_argument(
        "-r", "--reads", help="Reads SAM/BAM/FASTA/FASTQ file.", required=True
    )
    align_required_args.add_argument(
        "-s", "--segments", help="Segments FASTA/FASTQ file.", required=True
    )
    align_required_args.add_argument(
        "-b", "--boundaries", help="Sequence boundaries file.", required=True
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file in which to store alignment results.",
        default="extracted_bounded_sub_reads.fasta",
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
        "--prec_known",
        help=f"Probability of recombination for known segment alignment (default: {P_REC_KNOWN})",
        default=P_REC_KNOWN,
        type=float
    )

    parser.add_argument(
        "--prec_unknown",
        help=f"Probability of recombination for UNKNOWN segment alignment (default: {P_REC_UNKNOWN})",
        default=P_REC_UNKNOWN,
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
        help="Output file in which to store rejected reads.",
        default="extracted_bounded_sub_reads.rejected.fasta",
        required=False,
    )

    parser.add_argument(
        "--raw_marker_alignments",
        help="Output file in which to store raw alignments of markers for each section for each read.",
        default="extracted_bounded_sub_reads.raw_marker_alignments.txt",
        required=False,
    )

    parser.add_argument(
        "--initial_section_alignments",
        help="Output file in which to store initial section alignments per read.",
        default="extracted_bounded_sub_reads.initial_section_alignments.txt",
        required=False,
    )

    parser.add_argument(
        "--final_section_alignments",
        help="Output file in which to store filtered/final section alignments per read.",
        default="extracted_bounded_sub_reads.final_section_alignments.txt",
        required=False,
    )

    parser.add_argument(
        "-A",
        "--aligner",
        help="Use the given aligner.  [" + ", ".join([e.name for e in AlignmentAlgorithm]) + "]",
        type=str,
        required=False
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
    # args = parser.parse_args(args=raw_args)
    args = parser.parse_args()

    configure_logging(args)

    # Print logo:
    print_logo(args.aligner)

    # Log our command-line and log level so we can have it in the log file:
    LOGGER.info("Invoked by: %s", " ".join(raw_args))
    LOGGER.info("Complete runtime configuration settings:")
    for name, val in vars(args).items():
        LOGGER.info("    %s = %s", name, val)
    LOGGER.info("Log level set to: %s", logging.getLevelName(logging.getLogger().level))

    # Call our main method:
    extract_read_sections(args)

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == '__main__':
    main(sys.argv)
