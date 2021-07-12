import argparse
import gzip
from collections import OrderedDict
from math import ceil
import pysam
from construct import *

from multiprocessing.pool import ThreadPool
from functools import partial


def main():
    parser = argparse.ArgumentParser(description='Remove kinetics information', prog='remove_kinetics')
    parser.add_argument('-p', '--prefix', type=str, default="out", help="Output prefix")
    parser.add_argument('-x', '--exclude', type=str, default="ip,pw,fi,fp,ri,rp",
                        help='Comma-separated list of tags to exclude')
    parser.add_argument('bam', type=str, help="BAM")
    args = parser.parse_args()

    tags_to_exclude = args.exclude.split(",")

    # Silence message about the .bai file not being found.
    pysam.set_verbosity(0)

    bf = pysam.Samfile(args.bam, 'rb', check_sq=False)
    with pysam.Samfile(f'{args.prefix}.bam', 'wb', header=bf.header) as out:
        for read in bf:
            # Filter out unwanted tags
            if len(tags_to_exclude) > 0:
                filtered_tags = list(filter(lambda x: x[0] not in tags_to_exclude, read.get_tags()))
                read.set_tags(filtered_tags)

            out.write(read)


if __name__ == "__main__":
    main()
