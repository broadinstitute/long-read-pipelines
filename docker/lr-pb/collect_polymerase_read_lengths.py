import pysam
import argparse

import sys
import time
import logging
import re

################################################################################

LOGGER = logging.getLogger("collect_polymerase_read_lengths")

READ_NAME_RE = re.compile(r'.*?/(.*?)/(.*)')


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


def alignment_file_read_generator(open_alignment_file):
    for r in open_alignment_file.fetch(until_eof=True):
        yield r


def update_polymerase_read_count_with_next_read(read_gen, polymerase_read_length_map):
    try:
        cur_read = next(read_gen)
        zmw = int(READ_NAME_RE.match(cur_read.query_name).group(1))
        # NOTE: While positions start at 0, the end positions of the itervals are non-inclusive and indicate
        #       the length of the read, rather than the positions.
        max_pos = [int(b) for b in READ_NAME_RE.match(cur_read.query_name).group(2).split("_")][-1]
        try:
            if polymerase_read_length_map[zmw] < max_pos:
                polymerase_read_length_map[zmw] = max_pos
        except KeyError:
            polymerase_read_length_map[zmw] = max_pos
    except StopIteration:
        return True

    return False


def get_polymerase_read_lengths(args):
    """Get the lengths from the given polymerase reads and dump them to the specified output file."""

    with open(args.outfile, 'w') as out_file:

        out_file.write("\t".join(['zmw'] + ['polymerase_read_length']))
        out_file.write("\n")

        with pysam.AlignmentFile(args.subreads, 'rb', check_sq=False) as bam_file:

            # Account for optional scraps file:
            if args.scraps:
                scraps_file = pysam.AlignmentFile(args.scraps, 'rb', check_sq=False)

                scraps_generator = alignment_file_read_generator(scraps_file)
                all_scraps_processed = False
            else:
                all_scraps_processed = True

            try:
                total_num_reads = 0
                num_subreads = 0
                num_scraps = 0
                subreads_generator = alignment_file_read_generator(bam_file)

                polymerase_read_length_map = dict()

                all_subreads_processed = False

                while not all_subreads_processed or not all_scraps_processed:
                    # Get info for subreads:
                    if not all_subreads_processed:
                        all_subreads_processed = update_polymerase_read_count_with_next_read(subreads_generator,
                                                                                             polymerase_read_length_map)
                        num_subreads += 1
                        total_num_reads += 1

                    # Get info for scraps:
                    if not all_scraps_processed:
                        all_subreads_processed = update_polymerase_read_count_with_next_read(scraps_generator,
                                                                                             polymerase_read_length_map)
                        num_scraps += 1
                        total_num_reads += 1

                    if total_num_reads % 10000 == 0:
                        LOGGER.info(f"Processed {num_subreads} subreads\t{num_scraps} scraps")
            finally:
                if args.scraps:
                    scraps_file.close()

        LOGGER.info(f"Subreads processed: {num_subreads}")
        if args.scraps:
            LOGGER.info(f"Scraps processed: {num_scraps}")
        else:
            LOGGER.info(f"Total reads processed: {total_num_reads}")
        LOGGER.info(f"Writing results...")

        # Now write out the results:
        for k,v in polymerase_read_length_map.items():
            out_file.write(f"{k}\t{v}\n")

    LOGGER.info(f"Results written.")

################################################################################


def main(raw_args):
    base_outfile_name = "polymerase_read_lengths"

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        description="Ingests at least one file: "
                    "1) PacBio subreads file (SAM/BAM) containing raw subreads off the instrument.  "
                    "2) (optional) PacBio scraps file (SAM/BAM) containing raw scraps off the instrument.  "
                    "Uses the read names as indicators of polymerase read length (for speed).",
        usage="Count the lengths of the raw polymerase reads in each zmw.  "
              "If the scraps file is omitted, the resulting lengths will be ESTIMATES rather than true calculations.",
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-s", "--subreads", help="PacBio subreads SAM/BAM file.", required=True
    )
    required_args.add_argument(
        "-p", "--scraps", help="PacBio scraps SAM/BAM file.", required=False
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file in which to store results.",
        default=f"{base_outfile_name}.tsv",
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

    # Warn the user that without scraps, we can only do estimation:
    if not args.scraps:
        LOGGER.warning(f"No scraps file provided.  Counts produced will be APPROXIMATE!")

    # Call our main method:
    get_polymerase_read_lengths(args)

    LOGGER.info(f"Results written to: {args.outfile}")

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == '__main__':
    main(sys.argv)
