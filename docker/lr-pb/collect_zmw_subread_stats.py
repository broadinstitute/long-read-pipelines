import pysam
import argparse

import sys
import time
import logging
import re
import statistics

################################################################################

LOGGER = logging.getLogger("collect_zmw_subread_stats")


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


def get_stats(args):
    """Get the stats from the given polymerase reads and dump them to the specified output file(s)."""

    stat_names = ['num_subreads', 'min_length', 'max_length', 'range', 'mean', 'median', 'stdev']

    with open(args.outfile, 'w') as out_file:

        out_file.write("\t".join(['zmw'] + stat_names + ['raw_subread_lengths']))
        out_file.write("\n")

        with pysam.AlignmentFile(args.bam, 'rb', check_sq=False) as bam_file:

            zmw_re = re.compile(r'.*?/(.*?)/.*')

            prev_zmw = None
            zmw_read_lengths = []
            subread_seqs = []
            subread_names = []

            num_reads = 0

            for r in bam_file.fetch(until_eof=True):

                zmw = int(zmw_re.match(r.query_name).group(1))
                if prev_zmw and zmw != prev_zmw:

                    stats = dict()
                    stats['num_subreads'] = len(zmw_read_lengths)
                    stats["min_length"] = min(zmw_read_lengths)
                    stats["max_length"] = max(zmw_read_lengths)
                    stats["range"] = stats["max_length"] - stats["min_length"]
                    stats["mean"] = statistics.mean(zmw_read_lengths)
                    stats["median"] = statistics.median(zmw_read_lengths)

                    if len(zmw_read_lengths) > 1:
                        stats["stdev"] = statistics.stdev(zmw_read_lengths)
                    else:
                        stats["stdev"] = 0

                    out_file.write(f"{zmw}")
                    for k in stat_names:
                        if type(stats[k]) == int:
                            out_file.write(f"\t{stats[k]}")
                        else:
                            out_file.write(f"\t{stats[k]:.2f}")
                    out_file.write(f"\t{','.join(str(x) for x in zmw_read_lengths)}")
                    out_file.write("\n")

                    zmw_read_lengths = []
                    subread_seqs = []
                    subread_names = []

                zmw_read_lengths.append(len(r.query_sequence))
                subread_seqs.append(r.query_sequence)
                subread_names.append(r.query_name)
                prev_zmw = zmw

                num_reads += 1

                if num_reads % 5000 == 0:
                    LOGGER.info(f"Processed {num_reads} reads.")

        LOGGER.info(f"Total reads processed: {num_reads}")

################################################################################


def main(raw_args):

    base_outfile_name = "zmw_subread_stats"

    # Get our start time:
    overall_start = time.time()

    parser = argparse.ArgumentParser(
        description="Ingests one file: "
                    "PacBio subreads file (SAM/BAM) containing raw subreads off the instrument.",
        usage="Get the stats from the given polymerase reads and dump them to the specified output file(s).",
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-b", "--bam", help="PacBio subreads SAM/BAM file.", required=True
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

    # Call our main method:
    get_stats(args)

    LOGGER.info(f"Results written to: {args.outfile}")

    overall_end = time.time()
    LOGGER.info("Elapsed time: %f", overall_end - overall_start)


################################################################################


if __name__ == '__main__':
    main(sys.argv)
