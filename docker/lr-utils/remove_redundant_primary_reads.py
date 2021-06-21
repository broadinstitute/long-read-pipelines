import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description='Remove redundant reads from BAM file', prog='remove_redundant_reads')
    parser.add_argument('-p', '--prefix', type=str, default="shard", help="Output prefix")
    parser.add_argument('bam', type=str, help="BAM")
    args = parser.parse_args()

    # Silence message about the .bai file not being found.
    pysam.set_verbosity(0)

    seen_read_names = set()
    num_reads, num_primary_reads, num_redundant_reads = 0, 0, 0

    bf = pysam.Samfile(args.bam, 'rb', check_sq=False)
    with pysam.Samfile(f'{args.prefix}.bam', 'wb', header=bf.header) as out:
        for read in bf:
            num_reads += 1
            if read.is_secondary or read.is_supplementary:
                out.write(read)
            else:
                num_primary_reads += 1
                if read.query_name not in seen_read_names:
                    out.write(read)
                else:
                    num_redundant_reads += 1

                seen_read_names.add(read.query_name)

    print(f'num_reads\t{num_reads}')
    print(f'num_primary_reads\t{num_primary_reads}')
    print(f'num_redundant_reads\t{num_redundant_reads}')


if __name__ == "__main__":
    main()
