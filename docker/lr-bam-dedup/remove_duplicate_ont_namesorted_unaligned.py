import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description='Remove redundant reads from renamed-sorted ONT BAM file',
                                     prog='remove_redundant_reads')
    parser.add_argument('-p', '--prefix', type=str, default="shard", help="Output prefix")
    parser.add_argument('-q', '--qnames', type=str, help="Read names of duplicate records")

    parser.add_argument('bam', type=str, help="BAM")
    args = parser.parse_args()

    # Silence message about the .bai file not being found.
    pysam.set_verbosity(0)
    bf = pysam.Samfile(args.bam, 'rb', check_sq=False)

    num_records, num_dropped_records = 0, 0
    duplicate_record_names = list()

    with pysam.Samfile(f'{args.prefix}.bam', 'wb', header=bf.header) as out:

        # we rely on the observation that for queryname sorted, unaligned BAM,
        # if two neighboring records have the same query name, then they must be duplicate of each other
        current_qm = ''

        for read in bf:
            num_records += 1

            n = read.query_name
            if n == current_qm:
                duplicate_record_names.append(n)
                num_dropped_records += 1
            else:
                current_qm = n
                out.write(read)

    print(f'num_records: {num_records}')
    print(f'num_dropped_records: {num_dropped_records}')
    print(f'num_kept_alignments: {num_records - num_dropped_records}')

    with open(args.qnames, 'w') as outf:
        for qn in duplicate_record_names:
            outf.write(f'{qn}\n')


if __name__ == "__main__":
    main()
