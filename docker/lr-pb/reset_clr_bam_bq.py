import argparse
import pysam
import copy


def main():
    parser = argparse.ArgumentParser(description='Reset base qualities of reads in the CLR bam to the requested Phred base quality', prog='reset_clr_bam_bq')
    parser.add_argument('-q', '--basequal', type=int, default=10, help="Desired Phred base quality")
    parser.add_argument('-p', '--prefix', type=str, default="barbequed", help="Shard filename prefix")
    parser.add_argument('bam', type=str, help="BAM")
    args = parser.parse_args()

    # Silence message about the .bai file not being found.
    pysam.set_verbosity(0)

    if args.basequal < 0 or (args.basequal > 60 and args.basequal != 255):
        raise ValueError(f"Requested BQ value {args.basequal} isn't valid.")

    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.query_qualities
    bq = str(chr(args.basequal)) # no add 33 because the link above says so
    print(f"Setting base qualities to ASCII {str(chr(args.basequal+33))}.")

    bf = pysam.Samfile(args.bam, 'rb', check_sq=False)
    with pysam.Samfile(f'{args.prefix}.bam', 'wb', header=bf.header) as out:
        for read in bf:
            sausage = copy.deepcopy(read)
            n = len(sausage.query_sequence)
            qual = [ord(b) for b in list(bq * n)]
            sausage.query_qualities = qual
            out.write(sausage)


if __name__ == "__main__":
    main()