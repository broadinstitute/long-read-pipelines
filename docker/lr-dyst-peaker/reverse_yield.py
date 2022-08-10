import numpy as np
import argparse


######################################################################
def main():
    parser = argparse.ArgumentParser(description='Given a text file holding read lengths of all sequences, print lengths at which a certain fraction of reads are shorter than. The fraction bins are 10% to 90% with 10% increments.',
                                     prog='reverse_yield')
    parser.add_argument('-i', '--input', type=str, help="a text file holding read lengths of all sequences")
    parser.add_argument('-o', '--output', type=str, help="path to file holding length-9 reverse yield array (flat file)")
    args = parser.parse_args()

    with open(args.input) as inf:
        read_lengths = [int(line.strip()) for line in inf]

    sorted_read_lengths = sorted(read_lengths)

    with open(args.output, 'w') as outf:
        for frac in np.arange(0.1, 1, 0.1):
            outf.write( f"{sorted_read_lengths[round(frac*len(sorted_read_lengths))]}\n" )


######################################################################
if __name__ == "__main__":
    main()