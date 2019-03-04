#!/usr/bin/env python
from argparse import (ArgumentParser, FileType)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Convert Tandem Repeat Finder (TRF) dat file to bed format with repeat units for microsatellite genotyping')
    parser.add_argument(
        '--dat', type=str, required=True,
        help='Input dat file produced by Tandem Repeat Finder (TRF) using the -d option')
    parser.add_argument(
        '--bed', type=str, required=True,
        help='Output bed file containing genomic locations and repeat units of microsatellites.')

    return parser.parse_args()

### Main
def main():
    # Parse command line arguments
    args = parse_args()
    datfile = args.dat
    bedfile = args.bed

    with open(bedfile, 'w') as bed:
        chrom = ""
        with open(datfile, 'r') as dat:
            for line in dat:
                splitline = line.split()
                if line.startswith("Sequence:"):
                    chrom = line.split()[1]
                else:
                    # Catch index errors when line is blank
                    try:
                        # Check if in header sequence (all non-header lines start with an int: start pos)
                        try:
                            int(splitline[0])
                        except ValueError:
                            continue
                        start = splitline[0]
                        end = splitline[1]
                        motif = splitline[13]
                        copynum = splitline[3]
                        bed.write('\t'.join([chrom,start,end,motif,copynum]) + '\n')
                    except IndexError:
                        pass

if __name__ == '__main__':
    main()
