import pysam
import argparse

import sys
import time
import logging
import re
import statistics

################################################################################

LOGGER = logging.getLogger("add_pac_bio_rg_to_bam")


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


################################################################################


def add_read_group_to_bam(args):
    """Get the stats from the given polymerase reads and dump them to the specified output file(s)."""

    with pysam.Samfile(args.bam, 'rb', check_sq=False) as input_bam:

        header_dict = input_bam.header.as_dict()

        chars_to_skip = 0
        if args.header_read_group_line.startswith("@RG"):
            # Skip the first 4 characters (this includes the tab after the RG tag).
            chars_to_skip = 4
        elif args.header_read_group_line.startswith("RG"):
            # Skip the first 3 characters (this includes the tab after the RG tag).
            chars_to_skip = 3

        # Add the new header line to the header dict:
        rg_header_line = args.header_read_group_line[chars_to_skip:]
        header_dict["RG"] = rg_header_line

        new_header = pysam.AlignmentHeader.from_dict(header_dict)

        # Create a new output file:
        status_interval = 5000
        with pysam.AlignmentFile(args.outfile, 'wb', header=new_header) as out:
            for i, read in enumerate(input_bam):
                read.set_tag("RG", args.rg, value_type="i")
                out.write(read)

                if i % status_interval == 0:
                    print(f"Reads Processed:\t\t{i}")

        LOGGER.info(f"Total reads processed: {i}")

################################################################################


def main(raw_args):

    base_outfile_name = "out"

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        description="Requires one file and two arguments: "
                    "Reads file (SAM/BAM) containing reads to which to add the given read group ID.",
        usage="Add a PacBio read group to a bam file so it can be ingested with PacBio tools.",
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-b", "--bam", help="SAM/BAM file.", required=True
    )
    required_args.add_argument(
        "-r", "--rg", help="Read group ID.", required=True
    )
    required_args.add_argument(
        "-h", "--header-read-group-line", help="Read group header line to add to the bam header.", required=True
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="Output bam file with new read group info.",
        default=f"{base_outfile_name}.bam",
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

    # Log our command-line and log level so we can have it in the log file:
    LOGGER.info("Invoked by: %s", " ".join(raw_args))
    LOGGER.info("Complete runtime configuration settings:")
    for name, val in vars(args).items():
        LOGGER.info("    %s = %s", name, val)
    LOGGER.info("Log level set to: %s", logging.getLevelName(logging.getLogger().level))

    # Call our main method:
    add_read_group_to_bam(args)

    LOGGER.info(f"Results written to: {args.outfile}")

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == '__main__':
    main(sys.argv)
