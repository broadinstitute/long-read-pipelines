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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

import pysam

from tesserae import Tesserae
from tesserae.sequence import Sequence
from tesserae.tesserae import TesseraeAlignmentResult
from tesserae.tesserae import DEFAULT_REC

################################################################################

LOGGER = logging.getLogger("cartographer")

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


class SmithWatermanAlignmentResult:
    """
    Class to hold results from a SmithWaterman alignment.
    """

    def __init__(self, tgt_name, query_alignment_string, target_alignment_string,
                 query_start, query_end, target_start, target_end, score, max_score, qual_pl):
        self.name = tgt_name
        self.query_alignment_string = query_alignment_string
        self.target_alignment_string = target_alignment_string
        self.query_start = query_start
        self.query_end = query_end
        self.target_start = target_start
        self.target_end = target_end
        self.score = score
        self.max_score = max_score
        self.qual_pl = qual_pl

    def __str__(self):
        return (f"SmithWatermanAlignmentResult({self.name}, q[{self.query_start}-{self.query_end}], "
                f"t[{self.target_start}-{self.target_end}], {self.score} / {self.max_score}, {self.qual_pl})")


VALID_DNA_SEQUENCE_PATTERN = re.compile(r"^[ATGCNatgcn]+$")
MOSAIC_ALIGNER_LINE_MATCH_PATTERN = re.compile(r"^\s*(.*?)\s\[\s*(\d+)-\s*(\d+)\]\t+(.*)\s*$")
SMITH_WATERMAN_ALIGN_STRING_PATTERN = re.compile(r"^\s*(.*?)\s+\[pos:\s+(\d+);\s+len:\s+(\d+)\]\s*$")

DEFAULT_OUTPUT_FILE_NAME = "tesserae_alignment.bam"

DEFAULT_PLOT_FOLDER = "cartographer_plots"

# Default our quality to zero for now:
DEFAULT_BASE_QUALITY_CHAR = "!"

# Default max PL:
MAX_ALIGNMENT_PL = 60

# Minimum allowed PL quality value:
MIN_QUAL = 0

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


def log_variable(var_name, level=logging.DEBUG):
    """
    Logs the given variable's name and contents a the given log level.
    :param var_name: String containing the name of the variable to log.
    :param level: Level at which to log the given variable.
    :return: None
    """
    prev_frame = inspect.currentframe().f_back

    try:
        LOGGER.log(level, "%s = %s", var_name, prev_frame.f_locals[var_name])
    except KeyError:
        LOGGER.log(level, "%s = %s", var_name, prev_frame.f_globals[var_name])


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
            LOGGER.error("Error: You must provide a valid alignment algorithm.  Options are: %s", ", ".join([e.name  for e in AlignmentAlgorithm]))
            sys.exit(1)

    LOGGER.info("====================================================================")
    LOGGER.info("    ____           _                              _")
    LOGGER.info("   / ___|__ _ _ __| |_ ___   __ _ _ __ __ _ _ __ | |__   ___ _ __")
    LOGGER.info("  | |   / _` | '__| __/ _ \\ / _` | '__/ _` | '_ \\| '_ \\ / _ \\ '__|")
    LOGGER.info("  | |__| (_| | |  | || (_) | (_| | | | (_| | |_) | | | |  __/ |")
    LOGGER.info("   \\____\\__,_|_|   \\__\\___/ \\__, |_|  \\__,_| .__/|_| |_|\\___|_|")
    LOGGER.info("                            |___/          |_|")
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

    no_alignment = DetailedAlignmentInfo(-1, -1, 0, (), MIN_QUAL)

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
    return DetailedAlignmentInfo(start_index, end_index, template_length, cigar, aligned_qual_pl)


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
        if q < MIN_QUAL:
            return MIN_QUAL
        else:
            return q


def create_alignment_graph(query_alignment_info, target_aligned_info, target_sequences, title_prefix="", suffix="",
                           out_dir=DEFAULT_PLOT_FOLDER):
    """
    Create a graph of the alignment in the current working directory.

    :param query_alignment_info: A TesseraeAlignmentResult object representing the query string.
    :param target_aligned_info: A list[TesseraeAlignmentResult] representing the alignments of all targets to the
    query.
    :param target_sequences: Ordered list of tesserae.Sequence objects containing all targets to be aligned.
    :param title_prefix: Prefix to add to the title of the graph.
    :param suffix: Suffix to add to the output file.
    :param out_dir: Output folder into which to place the plots.  Must already exist.
    :return: None
    """

    # Create an easy way to look up our targets:
    target_name_set = set(t.seq_name for t in target_aligned_info)

    # Set up plotting variables before starting our plot:
    query_y = 1 + len(target_name_set)

    y_tick_pos = np.zeros(1 + len(target_name_set))
    y_tick_labels = []
    y_tick_pos[0] = query_y
    y_tick_labels.append("Query")

    target_dict = dict()
    target_num = 0
    for t in target_sequences:
        if t.name in target_name_set:
            target_dict[t.name] = (t, target_num)
            y_tick_pos[target_num + 1] = query_y - target_num - 1
            target_num += 1
            y_tick_labels.append(t.name)

    cmap = matplotlib.cm.get_cmap('Spectral')
    color_normalizer = matplotlib.colors.Normalize(vmin=0, vmax=len(target_name_set))

    font_dict = {"family": "monospace"}
    font_size = 6
    rect_border_width_px = 2

    # Center the rectangles at whole numbers on Y axis:
    rect_height = .75

    q_rect_height = .25
    q_rect_y_offset = -q_rect_height / 2

    xmin = 0
    xmax = len(query_alignment_info.alignment_string)
    ymin = 0
    ymax = query_y + 1

    fig_aspect_ratio = (16, 9)
    dpi = 200

    fig, axes = plt.subplots(num=1, figsize=fig_aspect_ratio, dpi=dpi)

    outline_color = [0, 0, 0, 1]
    query_color = [.7] * 3
    color_inv = [1 - c for c in query_color]

    # Plot the query:
    axes.add_patch(
        plt.Rectangle(
            xy=(xmin, query_y + q_rect_y_offset),
            width=len(query_alignment_info.alignment_string),
            height=q_rect_height,
            fill=True,
            linestyle="-",
            linewidth=rect_border_width_px,
            edgecolor=color_inv,
            facecolor=query_color
        )
    )

    # Plot the targets:
    for i, target in enumerate(target_aligned_info):
        target_num = target_dict[target.seq_name][1]

        y_coord = query_y - target_num - 1

        x_start = get_start_index_from_alignment_start_string(target.alignment_string)
        x_end = len(target.alignment_string) - 1

        color = cmap(color_normalizer(target_num))

        tgt_width = target.target_end_index - target.target_start_index

        x_rect_offset = ((x_end - x_start) / 2) - (tgt_width / 2)

        # Adjust the color by the quality:
        if isinstance(target, ProcessedAlignmentResult):
            qual_pl = target.overall_quality
        else:
            # Must calculate quality:
            detailed_alignment = compute_detailed_alignment_info(
                query_alignment_info.alignment_string,
                target.alignment_string,
                target.target_start_index,
                target.target_end_index
            )
            qual_pl = detailed_alignment.qual_pl
        qual_pl = math.floor(qual_pl)
        if qual_pl == 0:
            alpha = .1
            color = [.8] * 3 + [alpha]
            border_color = [.8] * 3 + [alpha]
            q_whisker_color = [.8] * 3 + [alpha]
        else:

            # By default the qual_pl is a log scale so the shading will be logarithmic:
            alpha = qual_pl / MAX_ALIGNMENT_PL

            # Uncomment below to do linear shading:
            raw_qual = math.pow(10, qual_pl / -10)
            alpha = 1-raw_qual

            _log_var(qual_pl)
            _log_var(raw_qual)
            _log_var(qual_pl / MAX_ALIGNMENT_PL)
            _log_var(alpha)

            color = [c for i, c in enumerate(color) if i < 3] + [alpha]
            border_color = [c for i, c in enumerate(outline_color) if i < 3] + [alpha]
            q_whisker_color = [c for i, c in enumerate(query_color) if i < 3] + [alpha]

        # Plot the alignment of the target:
        plot_target_query_position_box(
            axes,
            x_start,
            x_end,
            x_start + x_rect_offset,
            x_start + x_rect_offset + tgt_width,
            y_coord,
            rect_height,
            color,
            border_color,
            whisker_color=q_whisker_color,
            dir_width=x_rect_offset/4
        )

        # Plot the alignment qualities:
        axes.text(
            x_start + x_rect_offset + tgt_width/2 - tgt_width/8,
            y_coord - rect_height/2 - rect_height/8,
            f"Q{int(math.floor(qual_pl))}",
            color="k",
            fontsize=font_size,
            fontdict=font_dict
        )

        # Plot the alignment qualities:
        axes.text(
            x_start + x_rect_offset,
            query_y + rect_height/2 + rect_height/8,
            f"Q{int(math.floor(qual_pl))}",
            color="k",
            fontsize=font_size,
            fontdict=font_dict
        )

        # Plot the alignment of this target on the query:
        plot_target_query_position_box(
            axes,
            x_start,
            x_end,
            x_start + x_rect_offset,
            x_start + x_rect_offset + tgt_width,
            query_y,
            rect_height,
            color,
            border_color,
            whisker_color=color,
            dir_width=x_rect_offset/4
        )

    # Fix the plot's metadata:
    plt.axis([xmin, xmax, ymin, ymax])
    plt.title(f"{title_prefix} Alignment of Segments to {query_alignment_info.seq_name}")
    plt.ylabel("Sequence Number")
    plt.xlabel("Query Position (base number)")

    # Set y tick labels:
    plt.yticks(ticks=y_tick_pos, labels=y_tick_labels)

    clean_seq_name = query_alignment_info.seq_name.replace(os.path.sep, "_")
    pdf = PdfPages(out_dir + os.path.sep + f"alignments_{clean_seq_name}{suffix}.pdf")
    pdf.savefig(fig)
    pdf.close()

    plt.close(fig)


def plot_target_query_position_box(axes, query_x_start, query_x_end, tgt_x_start, tgt_x_end, y_coord, box_height,
                                   box_color, outline_color, whisker_color=[.7] * 3, dir_width=50):
    """
    Plot a box-and-whisker like box on the Y coordinate specified by `query_y` that represents the alignment between a
    target and a region inside a query (as specified by the arguments).
    :param axes: pyplot.axes object on which to plot the box.
    :param query_x_start: Left-most position of the alignment on the query.
    :param query_x_end: Right-most position of the alignment on the query.
    :param tgt_x_start: Left-most position of the alignment on the target.
    :param tgt_x_end: Right-most position of the alignment on the target.
    :param y_coord: Y coordinate on which to plot the box.
    :param box_height: Height of the box and whisker caps.
    :param box_color: Color of the fill of the box and of the whisker caps.
    :param outline_color: Color of the outline of the box.
    :param whisker_color: Color of the whiskers.
    :param dir_width: Width of the direction indicators on top of the caps.
    :return: None
    """

    # Plotting order is important to prevent occlusions and needing to fuss with the z-/layer value:

    # Left cap:
    plt.plot(
        [query_x_start, query_x_start],
        [y_coord + box_height/2, y_coord - box_height/2],
        '-',
        color=whisker_color
    )

    # Left cap direction indicators:
    for height_offset in [box_height/2, -box_height/2]:
        plt.plot(
            [query_x_start, query_x_start + dir_width],
            [y_coord + height_offset, y_coord + height_offset],
            '-',
            color=whisker_color
        )

    # Left whisker:
    plt.plot(
        [query_x_start, tgt_x_start],
        [y_coord, y_coord],
        '-',
        color=whisker_color
    )

    # Right whisker:
    plt.plot(
        [tgt_x_end, query_x_end],
        [y_coord, y_coord],
        '-',
        color=whisker_color
    )

    # Right cap:
    plt.plot(
        [query_x_end, query_x_end],
        [y_coord + box_height/2, y_coord - box_height/2],
        '-',
        color=whisker_color
    )

    # Right cap direction indicators:
    for height_offset in [box_height/2, -box_height/2]:
        plt.plot(
            [query_x_end - dir_width, query_x_end],
            [y_coord + height_offset, y_coord + height_offset],
            '-',
            color=whisker_color
        )

    # Box:
    axes.add_patch(
        plt.Rectangle(
            xy=(tgt_x_start, y_coord - box_height/2),
            width=tgt_x_end - tgt_x_start,
            height=box_height,
            fill=True,
            linestyle="-",
            linewidth=1,
            edgecolor=outline_color,
            facecolor=box_color
        )
    )


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

        # LOGGER.debug("Seq: %s", target_name)
        detailed_alignment_info = compute_detailed_alignment_info(
            read_result[1], sequence, target_start_index, target_end_index
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
            LOGGER.debug("Target does not pass quality or min base threshold: %s", target_name)
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


def find_all_ordered_seqs_in_alignment_results(seq_order, processed_results):
    """
    Search the given processed results for the given sequence order.

    Will search for the reverse-complemented read information as well as the given (assumed forward) direction
    information as specified by the RC_READ_NAME_IDENTIFIER suffix.

    :param seq_order: A list of strings representing the expected order of named sequences.
    :param processed_results: A list of ProcessedAlignmentResult containing alignments to search.
    :return: A list of dictionaries with one dictionary per complete sequence match.
             Each dict maps indices in the sequence order to the index in processed results.
             Partial matches are not returned.
             If the sequence is not contained in `processed_results`, returns an empty list.
    """

    hits = []

    # We need to account for unknown / barcode sequences here:
    base_seq_index = next(
        i for i, name in enumerate(seq_order) if name != UNKNOWN_NAME_IDENTIFIER and name != BARCODE_NAME_IDENTIFIER
    )

    # Create an ordering for the reverse complement:
    rev_seq_order = seq_order[::-1]

    match = dict()
    seq_index = base_seq_index
    for i, result in enumerate(processed_results):

        if result.seq_name == seq_order[seq_index] or \
            (result.seq_name.endswith(RC_READ_NAME_IDENTIFIER) and
             result.seq_name[:-len(RC_READ_NAME_IDENTIFIER)] == rev_seq_order[seq_index]):

            # Track our position in the processed_results if this is a new match:
            match[seq_index] = i

            if seq_index == len(seq_order) - 1:
                # Our match is complete!
                # Make sure we save it and reset our state.
                hits.append(match)
                match = dict()
                seq_index = base_seq_index
            else:
                # Get our next known sequence index:
                seq_index += 1
                while (seq_order[seq_index] == UNKNOWN_NAME_IDENTIFIER) \
                        or (seq_order[seq_index] == BARCODE_NAME_IDENTIFIER):
                    seq_index += 1
        else:
            # No match.  Reset our counters.
            seq_index = base_seq_index
            match.clear()

    return hits


def extract_meta_sequences_from_result(
        special_seq_id, result_index_dict, sequence_order, ordered_segment_alignment_results, read_sequence
):
    """
    Extract the special sequences as identified by special_seq_id from the given
    ordered_segment_alignment_results.

    Assumes the special_seq_id occurs in the given sequence_order.

    :param special_seq_id: String containing the special sequence ID to search for.
    :param result_index_dict: Dictionary mapping the index in sequence_order to the corresponding sequence in
    ordered_segment_alignment_results.
    :param sequence_order: The order of sequences to expect in ordered_segment_alignment_results.
    :param ordered_segment_alignment_results: A list of ProcessedAlignmentResult representing sequence alignments.
    :param read_sequence: The sequence of the read from which to extract the special sequences.
    :return: A list of MetaSequenceInfo objects, one for each special sequence discovered.
    """

    special_seq_indices = [i for i, n in enumerate(sequence_order) if n == special_seq_id]

    _log_var(result_index_dict)
    LOGGER.debug("Special seq indices in sequence order: %s", special_seq_indices)

    special_seq_info_list = []

    for indx in special_seq_indices:

        # This indx itself can never be in the result_index_dict because it corresponds to a special sequence that
        # we do not know ahead of time.
        # What should be guaranteed is that either or both of indx-1 and indx+1 will be in the result_index_dict.
        # Therefore we just check the corresponding conditions based on these values, rather than the index itself.

        if indx-1 not in result_index_dict:
            unknown_read_start_index = 0
        else:
            # The barcode should start right after the last base in the section before the barcode:
            tmp_index = result_index_dict[indx - 1]
            unknown_read_start_index = ordered_segment_alignment_results[tmp_index].read_end_pos + 1

        if indx+1 not in result_index_dict:
            unknown_read_end_index = len(read_sequence) - 1
        else:
            tmp_index = result_index_dict[indx + 1]
            # Subtract 1 here so we have correct inclusive coordinates for the sequences:
            unknown_read_end_index = ordered_segment_alignment_results[tmp_index].read_start_pos - 1

        # if indx == 0:
        #     unknown_read_start_index = 0
        # else:
        #     # The barcode should start right after the last base in the section before the barcode:
        #     tmp_index = result_index_dict[indx - 1]
        #     unknown_read_start_index = ordered_segment_alignment_results[tmp_index].read_end_pos + 1
        #
        # if result_index_dict[indx] == len(ordered_segment_alignment_results) - 1:
        #     unknown_read_end_index = len(read_sequence) - 1
        # else:
        #     tmp_index = result_index_dict[indx + 1]
        #     # Subtract 1 here so we have correct inclusive coordinates for the sequences:
        #     unknown_read_end_index = ordered_segment_alignment_results[tmp_index].read_start_pos - 1

        # Make sure to add 1 to the end position when pulling the sequence out so we account for the indexing using
        # inclusive coordinates.
        special_seq_info_list.append(
            MetaSequenceInfo(
                UNKNOWN_NAME_IDENTIFIER, read_sequence[unknown_read_start_index:unknown_read_end_index+1],
                unknown_read_start_index, unknown_read_end_index
            )
        )

    return special_seq_info_list


def dump_seq_map(seq_map, name='Reads'):
    """Dumps the given sequence map to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug('%s:', name)
        for e in seq_map.items():
            LOGGER.debug('    %s -> %s', e[0], e[1])


def dump_sequence_order(seq_list):
    """Dumps the given sequence list to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug('Ordered Sequences:')
        for i, s in enumerate(seq_list):
            LOGGER.debug('    % 2d: %s', i+1, s)


def dump_results(results_tuple_list):
    """Dumps the results tuple to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug('Results:')
        for r in results_tuple_list:
            LOGGER.debug('    %s', str(r))


def assert_valid_sequence_order(seq_list, known_segment_names_to_seq_dict):
    """
    Validate the given seq_list contains only known sequences in the given known_segment_names_to_seq_dict.

    To be valid, a barcode sequence can only occur once.

    :param seq_list: List of known sequences relating their expected relative order in each read.
    :param known_segment_names_to_seq_dict: dictionary of sequence name -> sequence upon which to validate seq_list.
    :return: None.  Raises a ValueError if the given seq_list is not valid.
    """

    has_seen_barcoode = False
    for i, s in enumerate(seq_list):
        if s not in known_segment_names_to_seq_dict and s != UNKNOWN_NAME_IDENTIFIER and s != BARCODE_NAME_IDENTIFIER:
            raise ValueError(
                f"Sequence {i+1} is not in the known sequence list: {s}.  "
                "Sequence names must either correspond to a sequence in the "
                f"given Segments file or be \"{UNKNOWN_NAME_IDENTIFIER}\" or \"{BARCODE_NAME_IDENTIFIER}\"",
            )
        elif s == BARCODE_NAME_IDENTIFIER:
            if has_seen_barcoode:
                raise ValueError(
                    "Given sequence order contains more than one barcode.  Sequence must contain at most one barcode."
                )
            else:
                has_seen_barcoode = True


def create_ordered_sequences(known_segments_to_seq_dict, sequence_order_list):
    """
    Create tesserae Sequence objects to align known ordered segments with tesserae.

    Adds both the given direction and the reverse-complemented direction so that the sequences can be aligned in either
    direction.  This is significantly faster than running the whole alignment twice.

    :param known_segments_to_seq_dict: dict of sequence_name -> sequence_bases
    :param sequence_order_list: list of sequence_name in the order in which they should appear in a read.
    :return: list(tesserae.Sequence) of all known sequences that are not UNKNOWN / BARCODE.
    """

    ordered_sequences = []

    for seq_name in sequence_order_list:
        if seq_name != UNKNOWN_NAME_IDENTIFIER and seq_name != BARCODE_NAME_IDENTIFIER:
            ordered_sequences.append(Sequence(seq_name, known_segments_to_seq_dict[seq_name]))
            ordered_sequences.append(
                Sequence(seq_name + RC_READ_NAME_IDENTIFIER, reverse_complement(known_segments_to_seq_dict[seq_name]))
            )

    return ordered_sequences


def create_unknown_sequence_info(sequence_order, known_segment_names_to_seq_dict):
    """
    Create tesserae Sequence objects for all sequences in the given dictionary that do not occur in the given
    sequence_order.

    This method in essence creates a list of sequences that could possibly be the UNKNOWN sections.

    :param sequence_order: List of sequence names in the order in which they should appear in a read.
    :param known_segment_names_to_seq_dict: dict of name -> sequence of all possible sequences to expect in a read.
    :return: list(tesserae.Sequence) of all known sequences that are not in the given sequence_order
    """

    sequences = []

    for name, seq in known_segment_names_to_seq_dict.items():
        if name not in sequence_order and name != UNKNOWN_NAME_IDENTIFIER and name != BARCODE_NAME_IDENTIFIER:
            sequences.append(Sequence(name, seq))

    return sequences


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


def ingest_barcode_file(barcode_file_path):
    """
    Read the contents of the given barcode file into a set of barcodes.

    Makes every effort to read the given file (stripping whitespace, removing quotes, ignoring any header lines).
    Performs basic validation on the barcodes as they are ingested and self-destructs if a barcode is invalid.

    :param barcode_file_path: Path to the barcode file.
    :return: set of strings - one for each barcode in the given barcode file.
    """

    barcode_set = set()

    with open(barcode_file_path, 'r') as f:
        is_not_first_line = False
        num_barcodes = 0
        for i, line in enumerate(f.readlines()):

            barcode = line.strip().replace("\"", "")

            if VALID_DNA_SEQUENCE_PATTERN.match(barcode):
                barcode_set.add(barcode)
                num_barcodes += 1
            elif is_not_first_line:
                LOGGER.error("Barcode %d is invalid: %s", i, barcode)
                LOGGER.error("CANNOT CONTINUE WITH INVALID BARCODES!  PLEASE FIX BARCODE FILE!")
                sys.exit(1)

            is_not_first_line = True

        LOGGER.debug("Read in %d barcodes.", num_barcodes)

    return barcode_set


def ingest_sequence_layout_order(file_path):
    """Ingest the contents of a layout file.

    A Layout/order file is a plain text file with one sequence name per line.
    These sequence names should correspond to the names of sequences in a FASTX
    file that is given as an argument to Cartographer.

    Names will have all whitespace at the start and end of the line stripped off.

    The validity of sequence names is not checked here."""

    seq_order = []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            seq_order.append(line.strip())

    return seq_order


def dump_results_to_log(results_tuple_list):
    """Dump the results_tuple_list to the log as a DEBUG message."""
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Results:")
        for result in results_tuple_list:
            LOGGER.debug("    %s", str(result))


################################################################################


def cartographize(args):
    """Main CLI call for the CARTographer tool."""

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

    # Create our plot folder if we must:
    if args.plot:
        LOGGER.info("Plots are enabled.  Will plot results for each alignment.")
        if os.path.exists(args.plot_folder):
            if os.path.isdir(args.plot_folder):
                LOGGER.info("Plot folder already exists.  Will overwrite plots in plot folder: %s", args.plot_folder)
            else:
                LOGGER.error("CANNOT CONTINUE.  A file exists in the plot folder location: %s", args.plot_folder)
                sys.exit(1)
        else:
            LOGGER.info("Creating plot folder (%s)...", args.plot_folder)
            os.makedirs(args.plot_folder)
    else:
        LOGGER.info("Plots are disabled.")

    LOGGER.info("Writing output to %s", args.outfile)
    if os.path.exists(args.outfile):
        LOGGER.warning("Outfile already exists.  Will overwrite: %s", args.outfile)

    structural_alignment_type = AlignmentAlgorithm.TESSERAE
    if args.aligner:
        try:
            structural_alignment_type = AlignmentAlgorithm[args.aligner]
        except KeyError:
            LOGGER.error("Error: You must provide a valid alignment algorithm.  Options are: %s",
                         ", ".join([e.name for e in AlignmentAlgorithm]))
            sys.exit(1)
    LOGGER.info("Structural Alignment Algorithm: %s", structural_alignment_type.name)

    LOGGER.info("Ingesting layout / order from %s ...", args.layout)
    sequence_order = ingest_sequence_layout_order(args.layout)
    LOGGER.info("Ingested %d ordered sequence names.", len(sequence_order))
    dump_sequence_order(sequence_order)

    # get number of unknown sequences:
    num_unknown_seqs = sequence_order.count(UNKNOWN_NAME_IDENTIFIER)
    LOGGER.info("Will search for %d %s sequences.", num_unknown_seqs, UNKNOWN_NAME_IDENTIFIER)

    # Get whether we're finding a barcode:
    do_search_for_barcode = BARCODE_NAME_IDENTIFIER in sequence_order
    if do_search_for_barcode:
        LOGGER.info("Will search for %s sequences.", BARCODE_NAME_IDENTIFIER)
    else:
        LOGGER.info("Will NOT search for %s sequences.", BARCODE_NAME_IDENTIFIER)

    if args.barcode_length:
        LOGGER.info("Expected barcode length: %d", args.barcode_length)

    barcode_set = None
    if args.barcode_file:
        if not args.barcode_length:
            LOGGER.error("Must always specify --barcode_length if specifying a barcode file.")
            sys.exit(1)
        LOGGER.info("Reading in barcode list from: %s", args.barcode_file)
        barcode_set = ingest_barcode_file(args.barcode_file)

    LOGGER.info("Ingesting known segments from %s ...", args.segments)
    known_segment_names_to_seq_dict, _ = ingest_fastx_file(args.segments)
    LOGGER.info("Ingested %d known segments.", len(known_segment_names_to_seq_dict))
    dump_seq_map(known_segment_names_to_seq_dict, "Known Segments")

    LOGGER.info("Validating given ordered sequences...")
    assert_valid_sequence_order(sequence_order, known_segment_names_to_seq_dict)
    LOGGER.debug("Sequence order is valid.")

    # Create an ordered dict of our ordered sequences:
    ordered_sequences = create_ordered_sequences(known_segment_names_to_seq_dict, sequence_order)

    # Create possible sequences in UNKNOWN sections:
    possible_unknown_sequences = create_unknown_sequence_info(sequence_order, known_segment_names_to_seq_dict)

    # A couple of spacing variables for nice looking logs:
    spacing_one = " " * 4
    spacing_two = " " * 8

    with open(args.outfile, 'w') as out_file:
        write_header_to_output_file(out_file, sequence_order, args.barcode_length, barcode_set)

        num_good_forward_alignments = 0
        num_good_reverse_complement_alignments = 0
        num_good_alignments = 0
        LOGGER.info("Processing reads...")

        with ReadFile(args.reads) as reads_file:
            num_reads = 0
            for read_num, read_sequence in enumerate(reads_file.get_reads()):
                num_reads += 1

                LOGGER.debug("%sProcessing read %d: %s (len=%d)",
                             spacing_one, read_num, read_sequence.name, len(read_sequence.sequence))

                index_dict_list = []

                LOGGER.info("%sInitial alignment of known segments with Tesserae...", spacing_one)
                ordered_segment_alignment_results, query_result = align_sequences(
                    read_sequence, ordered_sequences, args.minqual, args.minbases, p_rec=P_REC_KNOWN,
                    create_plots=args.plot, plot_folder=args.plot_folder, extra_plot_suffix="structural",
                    alignment_type=structural_alignment_type, threads=num_threads
                )

                if len(ordered_segment_alignment_results) != 0:
                    LOGGER.debug("%sDetermining if expected sequence appears in processed alignment results...",
                                 spacing_one)
                    index_dict_list = find_all_ordered_seqs_in_alignment_results(
                        sequence_order, ordered_segment_alignment_results
                    )

                LOGGER.info("%sExpected ordered sequence occurs %d time(s).", spacing_one, len(index_dict_list))

                if len(index_dict_list) == 0:
                    LOGGER.info("%sThis read has no complete matches for the given known ordered sequences.",
                                spacing_one)
                else:
                    LOGGER.info(
                        "%sAnalyzing read sub-sequences for %ss and %ss", spacing_one, UNKNOWN_NAME_IDENTIFIER,
                        BARCODE_NAME_IDENTIFIER
                    )

                    # Create a map from alignment string position to original read position:
                    read_alignment_string_original_position_map = create_alignment_to_base_map(
                        query_result.alignment_string
                    )

                    LOGGER.debug("Aligned/Original Pos Map: %s", read_alignment_string_original_position_map)

                    for result_index_dict in index_dict_list:

                        LOGGER.info(
                            "%sExtracting meta sequences from match at indices: %s", spacing_two,
                            sorted([i for i in result_index_dict.values()])
                        )

                        barcode = None
                        if do_search_for_barcode:
                            # At this point, we have enough information to extract the barcode if we must look for it:
                            barcode = extract_barcode(
                                ordered_segment_alignment_results,
                                query_result,
                                result_index_dict,
                                sequence_order,
                            )

                        unknown_seqs = []
                        if num_unknown_seqs != 0:

                            # Get our unknown sequences:
                            unknown_sequence_info = extract_unknown_sequences(
                                ordered_segment_alignment_results, query_result, result_index_dict, sequence_order
                            )

                            LOGGER.info(
                                "%sAligning %d remaining known sequences to %ss...", spacing_two,
                                len(possible_unknown_sequences), UNKNOWN_NAME_IDENTIFIER
                            )

                            for unknown_seq_num, unknown_seq_meta_alignment in enumerate(unknown_sequence_info):

                                # Clean any "insertions":
                                unknown_sequence = Sequence(
                                    f"{UNKNOWN_NAME_IDENTIFIER}{unknown_seq_num}",
                                    unknown_seq_meta_alignment.raw_read_alignment_string.replace("-", "")
                                )

                                LOGGER.debug(
                                    "%s%s %d [%d-%d]: %s", spacing_two, UNKNOWN_NAME_IDENTIFIER, unknown_seq_num,
                                    unknown_seq_meta_alignment.alignment_start_index,
                                    unknown_seq_meta_alignment.alignment_end_index,
                                    unknown_sequence
                                )

                                best_alignment = get_best_smith_waterman_alignment(
                                    unknown_sequence, possible_unknown_sequences
                                )

                                # Must change coordinates of the Smith-Waterman alignment to coorespond to original read
                                # and not the subsection of the read against which it was created:
                                best_alignment.query_start += unknown_seq_meta_alignment.alignment_start_index
                                best_alignment.query_end += unknown_seq_meta_alignment.alignment_start_index

                                LOGGER.info("%s%d = %s", UNKNOWN_NAME_IDENTIFIER, unknown_seq_num, best_alignment)
                                unknown_seqs.append(best_alignment)

                        known_sequence_alignemnts = \
                            [ordered_segment_alignment_results[j] for i, j in result_index_dict.items()]

                        # Note that regardless of the alignment direction we give the FORWARD sequence here:
                        write_entry_to_output_file(out_file, read_sequence.name, read_sequence.sequence,
                                                   is_forward_direction, read_num, known_sequence_alignemnts,
                                                   barcode, unknown_seqs, read_alignment_string_original_position_map,
                                                   args.barcode_length, barcode_set)

                        # Track how many good alignments we have:
                        num_good_alignments += 1

            LOGGER.info("Processed %d reads.", num_reads)
            LOGGER.info("# Forward segments with good structural alignments: %d", num_good_forward_alignments)
            LOGGER.info("# RC segments with good structural alignments: %d", num_good_reverse_complement_alignments)
            LOGGER.info("Total # reads with good structural alignments: %d", num_good_alignments)


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


def write_header_to_output_file(out_file, sequence_order, expected_barcode_length=None, barcode_set=None):
    """
    Write the header to the output file.  Assumes output file is in TSV format.
    :param out_file: File object to which to write.
    :param sequence_order: List containing sequence names and the order in which to find them in each read.
    :param expected_barcode_length: Expected length of a barcode.
    :param barcode_set: The set of known barcodes.
    :return: None
    """

    out_file.write("Read Name")
    out_file.write("\tRead Number")

    coord_specifier_string = "0-based, inclusive"

    for s in sequence_order:
        if s != UNKNOWN_NAME_IDENTIFIER and s != BARCODE_NAME_IDENTIFIER:
            out_file.write(f"\t{s} Start Pos ({coord_specifier_string})"
                           f"\t{s} End Pos ({coord_specifier_string})"
                           f"\t{s} Cigar"
                           f"\t{s} Map Qual")

    if BARCODE_NAME_IDENTIFIER in sequence_order:
        out_file.write(f"\tRaw {BARCODE_NAME_IDENTIFIER}"
                       f"\t{BARCODE_NAME_IDENTIFIER} Start Pos ({coord_specifier_string})"
                       f"\t{BARCODE_NAME_IDENTIFIER} End Pos ({coord_specifier_string})"
                       f"\t{BARCODE_NAME_IDENTIFIER} Length")
        if expected_barcode_length:
            out_file.write(f"\t{BARCODE_NAME_IDENTIFIER} With Context Bases")
        if barcode_set:
            out_file.write(f"\tOther Possible {BARCODE_NAME_IDENTIFIER}s")

    unk_count = 1
    for s in sequence_order:
        if s == UNKNOWN_NAME_IDENTIFIER:
            out_file.write(f"\t{UNKNOWN_NAME_IDENTIFIER}_{unk_count}"
                           f"\t{UNKNOWN_NAME_IDENTIFIER}_{unk_count} Read Start Pos ({coord_specifier_string})"
                           f"\t{UNKNOWN_NAME_IDENTIFIER}_{unk_count} Read End Pos ({coord_specifier_string})"
                           f"\t{UNKNOWN_NAME_IDENTIFIER}_{unk_count} Alignment Score"
                           f"\t{UNKNOWN_NAME_IDENTIFIER}_{unk_count} Max Alignment Score"
                           f"\t{UNKNOWN_NAME_IDENTIFIER}_{unk_count} Map Qual")
            unk_count += 1

    out_file.write("\tRead Sequence")

    out_file.write("\n")


def coord_transform(x, read_length, is_forward_direction=True):
    """
    Transforms the given coordinate, x, into forward direction coordinates if it is based in a reverse complemented
    read.

    Assumes x is 0-based.

    :param x: Coordinate to transform (in 0-based coordinates).
    :param read_length: Length of the read to which x refers.
    :param is_forward_direction: True iff the given coord, x, is from a reverse-complemented read.
    :return: The forward direction 0-based coordinate representation of x.
    """

    if is_forward_direction:
        return x

    return read_length - 1 - x


def write_entry_to_output_file(
        out_file,
        read_name,
        read_sequence,
        read_num,
        known_sequence_alignemnts,
        barcode,
        unknown_seqs,
        alignment_read_pos_map,
        expected_barcode_length=None,
        barcode_set=None
):
    """
    Write the alignment entry information to the given output file.

    :param out_file: Open File object to which to write.
    :param read_name: Name of the read to which the sequences were aligned.
    :param read_sequence: Sequence of the read to be aligned.
    :param read_num: Number of the read in the given FASTA file (0-based).
    :param known_sequence_alignemnts: List of ProcessedAlignmentResult, one for each known sequence in order.
    :param barcode: MetaSequenceInfo object containing alignment information about the barcode sequence (may be None).
    :param unknown_seqs: List of SmithWatermanAlignmentResult objects representing alignments for unknown sequences.
    :param alignment_read_pos_map: Dictionary mapping positions in the given alignment string back to the original
    string.  This is required to decode all positions from alignment string positions to original read positions.
    :param expected_barcode_length: The expected length of the barcode.
    :param barcode_set: The set of known barcodes.
    :return: None
    """

    out_file.write(read_name)
    out_file.write(f"\t{read_num + 1}")

    # TODO: Remove this hard-coded statement and replace on a section-by-section basis.
    is_forward_direction = True

    for a in known_sequence_alignemnts:
        start_coord, end_coord = get_well_formed_start_and_end_coords(
            alignment_read_pos_map[a.read_start_pos],
            alignment_read_pos_map[a.read_end_pos],
            len(read_sequence),
            is_forward_direction
        )

        out_file.write(f"\t{start_coord}"
                       f"\t{end_coord}"
                       f"\t{cigar_tuple_to_string(a.cigar)}"
                       f"\t{a.overall_quality}")

    if barcode:
        barcode_seq = barcode.raw_read_alignment_string.replace("-", "")
        barcode_length = len(barcode_seq)

        try:
            start_coord, end_coord = get_well_formed_start_and_end_coords(
                alignment_read_pos_map[barcode.alignment_start_index],
                alignment_read_pos_map[barcode.alignment_end_index],
                len(read_sequence),
                is_forward_direction
            )

            out_file.write(f"\t{barcode_seq}"
                           f"\t{start_coord}"
                           f"\t{end_coord}"
                           f"\t{barcode_length}")
        except KeyError:
            LOGGER.debug("Pos map: %s", alignment_read_pos_map)
            _log_var(barcode.alignment_start_index)
            _log_var(barcode.alignment_end_index)
            sys.exit(1)

        if expected_barcode_length:
            # If our barcode is the expected length, we don't need to do anything, but if it is less than that length,
            # we must include extra context around it:
            if barcode_length < expected_barcode_length:
                context = expected_barcode_length - barcode_length

                start_coord = (start_coord - context) if (start_coord - context) >= 0 else 0
                end_coord = (end_coord + 1 + context) \
                    if (end_coord + 1 + context) < len(read_sequence) \
                    else len(read_sequence) - 1

                fuzzy_barcode = read_sequence[start_coord:end_coord]
                if not is_forward_direction:
                    fuzzy_barcode = reverse_complement(fuzzy_barcode)

                out_file.write(f"\t{fuzzy_barcode}")
            else:
                out_file.write("\t-")

        if barcode_set:
            # OK, we must determine if we have a good barcode.
            # If so, we can output nothing here.
            # If not, we must give a list of likely barcodes here:
            if barcode_seq in barcode_set:
                out_file.write(f"\t{barcode_seq}")
            elif barcode_length < expected_barcode_length:
                out_file.write("\tMANY (Barcode length less than expected)")
            elif barcode_length > expected_barcode_length:
                possible_barcodes = []
                for i in range(barcode_length - expected_barcode_length + 1):
                    b = barcode_seq[i:i+expected_barcode_length]
                    if b in barcode_set:
                        possible_barcodes.append(b)

                if len(possible_barcodes) == 0:
                    out_file.write("\tMANY (Barcode subsequences not in barcode list.)")
                else:
                    out_file.write("\t")
                    out_file.write(",".join(possible_barcodes))
            else:
                out_file.write("\tMANY (Barcode not in list)")

    for unk in unknown_seqs:
        start_coord, end_coord = get_well_formed_start_and_end_coords(
            alignment_read_pos_map[unk.query_start],
            alignment_read_pos_map[unk.query_end],
            len(read_sequence),
            is_forward_direction
        )

        out_file.write(f"\t{unk.name}"
                       f"\t{start_coord}"
                       f"\t{end_coord}"
                       f"\t{unk.score}"
                       f"\t{unk.max_score}"
                       f"\t{unk.qual_pl}")

    out_file.write(f"\t{read_sequence}")
    out_file.write("\n")


def get_well_formed_start_and_end_coords(start_pos, end_pos, read_length, is_forward_direction):
    """
    Convert the given start and end positions to be correct for the strandedness of the read and return them in
    correct logical order (start < end).

    :param start_pos: Raw start position in the read.
    :param end_pos: Raw end Position in the read.
    :param read_length: Length of the read to which x refers.
    :param is_forward_direction: True iff the given coord, x, is from a reverse-complemented read.
    :return: A tuple containing a start and end position (0-based, inclusive) that are in order and are correct for the
    strandedness of the given read.
    """

    start_coord = coord_transform(start_pos, read_length, is_forward_direction)
    end_coord = coord_transform(end_pos, read_length, is_forward_direction)

    if start_coord > end_coord:
        tmp = end_coord
        end_coord = start_coord
        start_coord = tmp

    return start_coord, end_coord


def extract_unknown_sequences(ordered_segment_alignment_results, query_result, result_index_dict, sequence_order):
    """
    Extract the unknown sequences from the given alignment information.
    :param ordered_segment_alignment_results:
    :param query_result:
    :param result_index_dict:
    :param sequence_order:
    :return:
    """

    LOGGER.info("        Extracting unknown sequence info...")
    unknown_sequence_info = extract_meta_sequences_from_result(
        UNKNOWN_NAME_IDENTIFIER,
        result_index_dict, sequence_order, ordered_segment_alignment_results,
        query_result.alignment_string
    )
    LOGGER.info(
        "        Extracted %d %s sequences.", len(unknown_sequence_info), UNKNOWN_NAME_IDENTIFIER
    )
    return unknown_sequence_info


def extract_barcode(ordered_segment_alignment_results, query_result, result_index_dict, sequence_order):
    """
    Extract the barcode from the given alignments.

    :param ordered_segment_alignment_results:
    :param query_result:
    :param result_index_dict:
    :param sequence_order:
    :return:
    """

    LOGGER.info("        Extracting barcode info...")

    # Because there is only one barcode we can unwrap this immediately into an object:
    barcode_info = extract_meta_sequences_from_result(
        BARCODE_NAME_IDENTIFIER,
        result_index_dict, sequence_order, ordered_segment_alignment_results,
        query_result.alignment_string
    )[0]

    LOGGER.info(
        "        Barcode @ read pos %d-%d (len=%d): %s",
        barcode_info.alignment_start_index, barcode_info.alignment_end_index,
        len(barcode_info.raw_read_alignment_string), barcode_info.raw_read_alignment_string
    )
    return barcode_info


def create_result_from_smith_waterman_stdout(stdout, tgt_name):
    """
    Creates a TesseraeAlignmentResult object representing the alignment from the given stdout from a call to the
    SmithWaterman aligner.
    :param stdout: Stdout object (as created by subprocess.run) from which to create the alignment.
    :param tgt_name: Name of the target sequence being aligned.
    :return: A SmithWatermanAlignmentResult object representing the alignment in the given stdout.
    """

    # A Smith-Waterman record looks like this:
    # == Alignment 0 lengths (285, 45):
    #
    # hit 0.0 score: 20
    #   CTCAGCGAAGGAAAAACCCGCAGGAGG  [pos: 1; len: 27]
    #   ||||*| ||* ||*|||  |*||||||
    #   CTCAAC-AAC-AACAAC--GGAGGAGG  [pos: 6; len: 23]
    #
    # ==

    lines = stdout.decode("utf-8").split("\n")

    # Get the lengths of the aligned strings:
    needle = "lengths ("
    q_len, t_len = lines[0][lines[0].find(needle) + len(needle):-1][:-1].strip().split(', ')
    q_len = int(q_len)
    t_len = int(t_len)

    _log_var(q_len)
    _log_var(t_len)

    max_score = min(q_len, t_len) * 2

    # Get the score
    needle = "score:"
    score = int(lines[2][lines[2].rfind(needle)+len(needle):].strip())

    # Get the query info:
    query_match = SMITH_WATERMAN_ALIGN_STRING_PATTERN.match(lines[3])
    if not query_match:
        LOGGER.error("Could not parse SmithWaterman results for query from line: %s", lines[3])
        sys.exit(1)
    query_align_string = query_match.group(1)
    query_start = int(query_match.group(2))
    query_end = query_start + int(query_match.group(3))

    # Get the target info:
    target_match = SMITH_WATERMAN_ALIGN_STRING_PATTERN.match(lines[5])
    if not target_match:
        LOGGER.error("Could not parse SmithWaterman results for target from line: %s", lines[5])
        sys.exit(1)
    target_align_string = target_match.group(1)
    target_start = int(target_match.group(2))
    target_end = target_start + int(target_match.group(3))

    # Get the number of errors in the sequence:
    l = lines[4].strip()
    num_errors = l.count("*") + l.count(" ")

    _log_var(num_errors)
    _log_var(len(query_align_string))
    _log_var(len(target_align_string))

    return SmithWatermanAlignmentResult(
        tgt_name,
        query_align_string,
        target_align_string,
        query_start,
        query_end,
        target_start,
        target_end,
        score,
        max_score,
        get_qual_pl(num_errors, len(query_align_string))
    )


def get_best_smith_waterman_alignment(read_sequence, target_sequences, log_spacing="    "):
    """
    Perform a Smith Waterman alignment on the given data and return back  the single best alignment.

    :param read_sequence: tesserae Sequence object against which to align all sequences in target_sequences
    :param target_sequences: ordered list of tessereae Sequence objects to align.
    :param log_spacing: Spacing to precede any log statements.
    :return: A list of SmithWatermanAlignmentResult objects.
    """

    best = None

    _, seq_file = tempfile.mkstemp()
    try:
        for t in target_sequences:
            smith_waterman_args = [
                "/seq-align/bin/smith_waterman",
                "--minscore", '0',
                "--maxhits", '1',
                "--pretty",
                read_sequence.sequence,
                t.sequence
            ]

            LOGGER.debug(
                "%sRunning Smith-Waterman on tmp file (%s): %s", log_spacing, seq_file, " ".join(smith_waterman_args)
            )

            completed_process = subprocess.run(smith_waterman_args, capture_output=True, check=True)

            if LOGGER.isEnabledFor(logging.DEBUG):
                LOGGER.debug("=" * 80)
                LOGGER.debug("Raw SmithWaterman Output:")
                LOGGER.debug("%s", completed_process.stdout.decode("utf-8"))
                LOGGER.debug("=" * 80)

            r = create_result_from_smith_waterman_stdout(completed_process.stdout, t.name)
            LOGGER.debug("Smith-Waterman Results: %s", r)
            LOGGER.debug("")

            # We use the scores here rather than the qualities because qual is capped and score is not.
            if not best:
                best = r
            elif (best.score/best.max_score) < (r.score/r.max_score):
                best = r
    finally:
        os.remove(seq_file)

    return best


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
        _ = subprocess.run(mosaic_aligner_args, capture_output=True, check=True)

        if LOGGER.isEnabledFor(logging.DEBUG):
            with open(out_file_name, 'rb') as f:
                LOGGER.debug("=" * 80)
                LOGGER.debug("Raw MosaicAligner Alignment:")
                for line in f.readlines():
                    LOGGER.debug("%s", line.decode("ascii").rstrip())
                LOGGER.debug("=" * 80)

        return get_raw_results_from_mosaic_aligner_file()

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
                    p_rec=DEFAULT_REC, alignment_type=AlignmentAlgorithm.TESSERAE, create_plots=False,
                    plot_folder=DEFAULT_PLOT_FOLDER, extra_plot_suffix="",
                    threads=1, log_spacing="    "):
    """
    Perform a Tesserae alignment on

    :param read_sequence: tesserae Sequence object against which to align all sequences in target_sequences
    :param target_sequences: ordered list of tessereae Sequence objects to align.
    :param minqual: Minimum quality for an alignment to be retained.
    :param minbases: Minimum number of bases for an alignment to be retained.
    :param p_rec: Prior probability of recombination event taking place.
    :param alignment_type: AlignmentAlgorithm object specifying which alignment algorithm to use.
    :param create_plots: If true, will create plots for the resulting alignment in the given plot_folder.
    :param plot_folder: Folder into which to place plots.  Must already exist.
    :param extra_plot_suffix: String to append to the end of the plot name.
    :param threads: Number of threads with which to run Tesserae alignment.
    :param log_spacing: Spacing to precede any log statements.
    :return: A tuple of results and the read sequence result (list(ProcessedAlignmentResult), TesseraeAlignmentResult)
    """

    start_time = time.time()
    if alignment_type == AlignmentAlgorithm.MOSAIC_ALIGNER:
        raw_results = create_raw_alignment_with_mosaic_aligner(read_sequence, target_sequences)
    elif alignment_type == AlignmentAlgorithm.SMITH_WATERMAN:
        LOGGER.error("SMITH WATERMAN not implemented yet for structural alignment.")
        sys.exit(1)
        # raw_results = create_raw_alignment_with_smith_waterman(read_sequence, target_sequences)
    else:
        raw_results = Tesserae(prho=p_rec, threads=threads).align(read_sequence, target_sequences)
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

        if create_plots:
            # Try to make a PDF:
            if alignment_type == AlignmentAlgorithm.MOSAIC_ALIGNER:
                decorator = "_mosaicAligner"
            elif alignment_type == AlignmentAlgorithm.SMITH_WATERMAN:
                decorator = "_smithWaterman"
            else:
                decorator = "_tesserae"
            create_alignment_graph(
                query_result, raw_results[1:],
                target_sequences,
                title_prefix="Raw " + extra_plot_suffix.capitalize(),
                suffix="_" + extra_plot_suffix + "_raw" + decorator,
                out_dir=plot_folder,
            )
            create_alignment_graph(
                query_result,
                alignment_results,
                target_sequences,
                title_prefix="Filtered " + extra_plot_suffix.capitalize(),
                suffix="_" + extra_plot_suffix + "_filtered" + decorator,
                out_dir=plot_folder,
            )

    return alignment_results, query_result


################################################################################


def main(raw_args):

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        description="Ingests three files: "
                    "one read file(FASTA/FASTQ/SAM/BAM) containing reads to be mapped, "
                    "one FASTA file containing known possible sequences that can occur in the read."
                    "and a file containing possible sequence order.  This sequence order file is"
                    "a plain text file with one sequence name per line.  The names should correspond"
                    "to the sequence names in the given sequence FASTA file.  If a sequence is not "
                    f"known it should be labeled {UNKNOWN_NAME_IDENTIFIER}.  If a sequence is a barcode it should be "
                    f"labeled {BARCODE_NAME_IDENTIFIER}.",
        usage="map out the order of sequences in a set of reads",
    )

    align_required_args = parser.add_argument_group("required arguments")
    align_required_args.add_argument(
        "-r", "--reads", help="Reads SAM/BAM/FASTA/FASTQ file.", required=True
    )
    align_required_args.add_argument(
        "-s", "--segments", help="Segments FASTA/FASTQ file.", required=True
    )
    align_required_args.add_argument(
        "-l", "--layout", help="Sequence layout / order file.", required=True
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file in which to store alignment results.",
        default="cartograph_output.tsv",
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
        "--barcode_length",
        help="Length of the barcode sequence",
        type=int,
        required=False
    )

    parser.add_argument(
        "--barcode_file",
        help="File containing a list of all possible barcodes.  If specified, must also specify --barcode_length.",
        type=str,
        required=False
    )

    parser.add_argument(
        "-A",
        "--aligner",
        help="Use the given aligner.  [" + ", ".join([e.name for e in AlignmentAlgorithm]) + "]",
        type=str,
        required=False
    )

    parser.add_argument(
        "--plot",
        help=f"If true, will create plots for each alignment in the plot folder.",
        action="store_true"
    )

    parser.add_argument(
        "--plot_folder",
        help=f"Folder in which to put plots (default: {DEFAULT_PLOT_FOLDER})",
        default=DEFAULT_PLOT_FOLDER,
        type=str
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
    cartographize(args)

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == '__main__':
    main(sys.argv)
