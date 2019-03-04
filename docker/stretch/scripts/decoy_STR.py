#!/usr/bin/env python
"""Generate a fasta file containing all possible STR repeat units.
This can be used as a decoy for STR reads to map to.
"""

import sys
from argparse import (ArgumentParser, FileType)
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import itertools

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = ArgumentParser(description='Generate a fasta file containing all ' \
        'possible STR repeat units.')
    parser.add_argument(
        '--output', type=str, required=False,
        help='Output filename for fasta file. Defualt: stdout.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '--length', type=int, required=False, default=500,
        help='Length of repetitive sequence in bp to simulate for each decoy locus.')
    group.add_argument(
        '--repeatunits', type=int, required=False,
        help='Number of repeat units to simulate for each decoy locus.')


    return parser.parse_args()

def circular_permuted(x):
    return([x[i:] + x[:i] for i in range(len(x))])

def self_and_rev_complement(in_dna):
    all_possible = [in_dna]

    # Get reverse complement
    dna = Seq(in_dna, generic_dna)
    rev_complement = str(dna.reverse_complement())
    all_possible.append(rev_complement)
    return(all_possible)

def normalise_str(in_dna):
    """Find all possible eqivalent STR sequences.
    And return the first alphabetically.

    For example, TA = AT. But would return AT.
    """
    all_possible = []
    # Circularly permute original sequence and reverse complement
    for seq in self_and_rev_complement(in_dna):
        for permuted_seq in circular_permuted(seq): # Switch to faster permutation (6)
            all_possible.append(permuted_seq)

    # Sort and take the first
    all_possible.sort()
    return(all_possible[0])

def main():
    # Parse command line arguments
    args = parse_args()
    if args.output:
        outstream = open(args.output, 'w')
    else:
        outstream = sys.stdout
    seqlength = args.length


    args.length

    nucleotides = 'ATCG'
    repeatlenrange = range(1,7) #1-6
    allRUs = []
    existingRUs = set()
    for i in repeatlenrange:
        newRUs = set([normalise_str(''.join(x)) for x in itertools.product(nucleotides, repeat=i)])
        existingRUs.update([''.join(itertools.repeat(j, i)) for j in allRUs for i in repeatlenrange])

        allRUs.extend([RU for RU in newRUs if RU not in existingRUs])

    # Write fasta file
    for RU in allRUs:
        RUlen = len(RU)
        if args.repeatunits:
            nrepeats = args.repeatunits
            sequence = ''.join(itertools.repeat(RU, nrepeats))
        else:
            nrepeats = int(seqlength / RUlen + RUlen)
            sequence = ''.join(itertools.repeat(RU, nrepeats))[:seqlength]

        record = SeqRecord(Seq('{}'.format(sequence),
                        generic_dna), id='STR-{}'.format(RU), description = '')
        SeqIO.write(record, outstream, "fasta")

if __name__ == '__main__':
    main()
