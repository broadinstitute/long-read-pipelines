import pysam
import copy
import argparse


def reverse_orig_seq_if_necessary(input_aln: pysam.AlignedSegment) -> pysam.AlignedSegment:
    if input_aln.is_mapped and input_aln.is_reverse:
        seq = input_aln.get_forward_sequence()
        bases_q = input_aln.get_forward_qualities()

        input_aln.query_sequence = seq
        input_aln.query_qualities = bases_q


def main():
    parser = argparse.ArgumentParser(description='Reset the fjsl',
                                     prog='reset_alignment')
    parser.add_argument('-i', '--input', type=str, required=True, help="input bam name")
    parser.add_argument('-o', '--output', type=str, required=True,  help="Output bam name")
    args = parser.parse_args()

    mode = 'rb' if args.input.endswith('.sam') else 'r'
    input_file = pysam.AlignmentFile(args.input, mode)
    output_file = pysam.AlignmentFile(args.output, 'wb', header=input_file.header)
    try:
        cnt = 0
        for orig_aln in input_file.fetch():
            reverse_orig_seq_if_necessary(orig_aln)
            orig_aln.flag = 4  # mark all as unmapped
            cnt += 1
            output_file.write(orig_aln)
        print(f"Total records processed: {cnt}.")
    finally:
        output_file.close()
        input_file.close()

if __name__ == "__main__":
    main()