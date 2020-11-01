import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description='Extract uncorrected reads from a subreads BAM',
                                     prog='extract_uncorrected_reads')
    parser.add_argument('-o', '--output', type=str, help="output BAM")
    parser.add_argument('subreads_bam', type=str, help="subreads BAM")
    parser.add_argument('consensus_bam', type=str, help="corrected BAM")
    args = parser.parse_args()

    # Populate a set of all ZMWs seen in the consensus BAM.  Note that in Python,
    # each `int` is actually 24 bytes.  Assuming the worst case scenario of a unique
    # ZMW for 8-million ZMWs in a SMRTCell 8M, this should be ~192 megabytes of storage.
    cf_zmws = set()
    cf = pysam.Samfile(args.consensus_bam, 'rb', check_sq=False)
    for cf_read in cf:
        cf_zmw = cf_read.get_tag("zm")
        cf_zmws.add(cf_zmw)

    # Now iterate through the subreads, and if a ZMW is not in the corrected list,
    # emit the read to the output file.
    num_uncorrected = 0
    sf = pysam.Samfile(args.subreads_bam, 'rb', check_sq=False)
    with pysam.Samfile(args.output, 'wb', header=sf.header) as out:
        for sf_read in sf:
            sf_zmw = sf_read.get_tag("zm")

            if sf_zmw not in cf_zmws:
                out.write(sf_read)
                num_uncorrected += 1

    print(f'Wrote {num_uncorrected} uncorrected reads to {args.output}.')


if __name__ == "__main__":
    main()
