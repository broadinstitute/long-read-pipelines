import numpy as np
import argparse
from scipy.stats import skew

######################################################################
def main():
    parser = argparse.ArgumentParser(description='Given a text file holding read lengths of all sequences, print G1 skewness.',
                                     prog='measure_g1_skew')
    parser.add_argument('-i', '--input', type=str, help="a text file holding read lengths of all sequences")
    parser.add_argument('-o', '--output', type=str, help="path to output file")
    args = parser.parse_args()

    with open(args.input) as inf:
        read_lengths = [int(line.strip()) for line in inf]

    with open(args.output, 'w') as outf:
        outf.write(f'{skew(read_lengths):.2f}\n')


######################################################################
if __name__ == "__main__":
    main()