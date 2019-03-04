#!/usr/bin/env python
"""Sort a fasta file alphabetically by chromosome name/sequence ID
"""

import argparse
import sys
from Bio import SeqIO

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Sort a fasta file alphabetically by chromosome name/sequence ID.')
    parser.add_argument(
        '--infile', type=str, required=True,
        help='Input fasta file.')
    parser.add_argument(
        '--outfile', type=str, required=True,
        help='Output fasta filename.')
    return parser.parse_args()

fastafile = "/group/bioi1/harrietd/git/STR-pipelines/simulate_reads/reference-data/hg19.STRdecoys.fasta"
outfasta = "hg19.STRdecoys.sorted.fasta"


def main():
    # Parse command line arguments
    args = parse_args()
    infile = args.infile
    outfile = args.outfile

    with open(outfile, "w") as outfasta:

        allfasta = SeqIO.parse(infile, "fasta")
        ids = sorted([rec.id for rec in allfasta])
        record_index = SeqIO.index(infile, "fasta")
        records = (record_index[id] for id in ids)
        SeqIO.write(records, outfile, "fasta")

if __name__ == '__main__':
    main()
