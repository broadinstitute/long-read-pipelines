#!/usr/bin/env python3
"""
Filter annotations in a GFF file that are predicted to be on the chromosome, using MOB-suite's
contig classifications
"""

import argparse
from pathlib import Path

import pandas


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.strip()
    )

    parser.add_argument(
        'gff', type=Path, help="Path to GFF file to filter"
    )
    parser.add_argument(
        'contig_report', type=Path, help="Path to MOB-suite contig report"
    )

    args = parser.parse_args()

    contig_report = pandas.read_csv(args.contig_report, sep='\t')

    to_keep = set(map(lambda v: v.split()[0], contig_report.query('molecule_type == "plasmid"')['contig_id']))

    in_fasta = False
    with args.gff.open() as f:
        for line in f:
            line = line.strip()
            if in_fasta:
                print(line)
            else:
                if line.startswith('#'):
                    if line == "##FASTA":
                        print(line)
                        in_fasta = True
                    elif line.startswith('##sequence-region'):
                        parts = line.split()
                        if parts[1] in to_keep:
                            print(line)
                    else:
                        print(line)
                else:
                    parts = line.split('\t')
                    if parts[0] in to_keep:
                        print(line)


if __name__ == '__main__':
    main()
