import pysam
import argparse

import sys
import time
import logging
import re

################################################################################

LOGGER = logging.getLogger("add_pac_bio_rg_to_bam")

HEADER_SPLITTER = re.compile(r'''[ \t]+''')

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

    hex_chars = [chr(i+ord('a')) for i in range(6)] + [chr(i+ord('A')) for i in range(6)]
    if any(c in args.rg for c in hex_chars):
        new_rg = int(args.rg, 16)
    else:
        new_rg = args.rg

    with pysam.AlignmentFile(args.bam, 'rb', check_sq=False) as input_bam:

        header_text = input_bam.text
        print('-' * 80)
        print(header_text)
        print('-' * 80)

        chars_to_skip = 0
        if args.header_read_group_line.startswith("@RG"):
            # Skip the first 4 characters (this includes the tab after the RG tag).
            chars_to_skip = 4
        elif args.header_read_group_line.startswith("RG"):
            # Skip the first 3 characters (this includes the tab after the RG tag).
            chars_to_skip = 3

        # grab the contents of the read group line:
        raw_rg_header_line = args.header_read_group_line[chars_to_skip:]

        # Do a little cleanup here to account for copy / paste and input errors:
        # Split the line and create entries in sub-dict for each field:
        rg_field_dict = dict()
        for field_value_string in HEADER_SPLITTER.split(raw_rg_header_line):
            if len(field_value_string) == 0:
                continue
            i = field_value_string.find(":")
            field_name = field_value_string[:i]
            if field_name == "ID":
                field_value = new_rg
            else:
                field_value = field_value_string[i+1:]
            rg_field_dict[field_name] = field_value

        clean_rg_header_line = "@RG\t" + "\t".join([f"{k}:{v}" for k, v in rg_field_dict.items()])

        new_header = pysam.AlignmentHeader.from_text(header_text + clean_rg_header_line)

        # Create a new output file:
        status_interval = 5000
        with pysam.AlignmentFile(args.outfile, 'wb', header=new_header) as out:
            print('-' * 80)
            print(out.text)
            print('-' * 80)
            for i, read in enumerate(input_bam):
                read.set_tag("RG", new_rg, value_type="i")
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
        "-l", "--header-read-group-line", help="Read group header line to add to the bam header.", required=True
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
