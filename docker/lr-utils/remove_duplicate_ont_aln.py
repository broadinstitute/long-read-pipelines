import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description='Remove redundant alignment records from ONT BAM file',
                                     prog='remove_redundant_reads')
    parser.add_argument('-p', '--prefix', type=str, default="shard", help="Output prefix")
    parser.add_argument('-a', '--annotations', type=str, help="Annotations on (potential) duplicate reads")

    parser.add_argument('bam', type=str, help="BAM")
    args = parser.parse_args()

    # create a dict of set's, a trick to avoid Hash collisions
    guilty_dict_per_chr = dict()
    with open(args.annotations) as f:
        for line in f:
            arr = line.strip().split('\t')
            name = arr[0]
            chrom = arr[2]
            guilty_dict_per_chr.setdefault(chrom, set())
            guilty_dict_per_chr[chrom].add(name)

    print("chromosomes on which there are duplicate records:")
    print(f"  {guilty_dict_per_chr}")

    # Silence message about the .bai file not being found.
    pysam.set_verbosity(0)

    num_alignments, num_dropped_alignments = 0, 0
    bf = pysam.Samfile(args.bam, 'rb', check_sq=False)
    with pysam.Samfile(f'{args.prefix}.bam', 'wb', header=bf.header) as out:
        # we rely on the observation that for coordinate sorted BAM,
        # duplicate records will appear in blocks, hence once we step off a position with duplicates, we start afresh
        current_position = -1
        current_signatures = set()
        for read in bf:
            num_alignments += 1

            chrom = read.reference_name
            n = read.query_name
            if chrom in guilty_dict_per_chr and n in guilty_dict_per_chr[chrom]:

                mq = read.mapping_quality
                sam_flag = read.flag
                pos = read.reference_start
                cigar = read.cigarstring
                signature = f"{n}-{chrom}-{pos}-{mq}-{sam_flag}-{cigar}"

                if current_position != pos:  # new position, let's write and reset
                    out.write(read)
                    current_position = pos
                    current_signatures = set()
                    current_signatures.add(signature)
                elif signature in current_signatures:  # You're wanted!
                    num_dropped_alignments += 1
                    pass
                else:  # you are a new group of duplicates that map to this location
                    out.write(read)
                    current_signatures.add(signature)
            else:
                out.write(read)

    print(f'num_alignments: {num_alignments}')
    print(f'num_dropped_alignments: {num_dropped_alignments}')
    print(f'num_kept_alignments: {num_alignments - num_dropped_alignments}')


if __name__ == "__main__":
    main()
