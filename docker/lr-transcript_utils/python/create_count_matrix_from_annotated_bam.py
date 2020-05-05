#!/usr/bin/env python

import os
import sys
import argparse

import pysam
import tqdm

gGENE_TX_TAG_NAME = "XG"
gCBC_TAG_NAME = "CB"
gUMI_TAG_NAME = "BX"


def main(input_bam, out_tsv_name):

    global gGENE_TX_TAG_NAME
    global gCBC_TAG_NAME
    global gUMI_TAG_NAME

    print("Verifying input file(s) exist...", file=sys.stderr)
    files_ok = True
    for f in [input_bam]:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}", file=sys.stderr)
            files_ok = False
    if not files_ok:
        sys.exit(1)
    print("Input files verified.", file=sys.stderr)

    reads_processed = 0

    print(f"Tallying gene/tx x CBC x UMI counts from {input_bam}", file=sys.stderr)
    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
            input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        file=sys.stderr
    ) as pbar:

        # Create a storage system for our count matrix:
        count_matrix = dict()

        has_cbc = True

        for read in bam_file:

            gene_tx = read.get_tag(gGENE_TX_TAG_NAME)
            umi = read.get_tag(gUMI_TAG_NAME)

            if has_cbc:
                try:
                    cbc = read.get_tag(gCBC_TAG_NAME)
                except KeyError:
                    has_cbc = False

            if not has_cbc:
                cbc = "SINGLE_CELL"

            if cbc in count_matrix:
                if gene_tx in count_matrix[cbc]:
                    if umi in count_matrix[cbc][gene_tx]:
                        count_matrix[cbc][gene_tx][umi] += 1
                    else:
                        count_matrix[cbc][gene_tx][umi] = 1
                else:
                    count_matrix[cbc][gene_tx] = dict()
                    count_matrix[cbc][gene_tx][umi] = 1
            else:
                count_matrix[cbc] = dict()
                count_matrix[cbc][gene_tx] = dict()
                count_matrix[cbc][gene_tx][umi] = 1

            # Update our counts:
            pbar.update(1)
            reads_processed += 1

    print("Done!", file=sys.stderr)
    print(f"Reads processed: {reads_processed}", file=sys.stderr)

    print(f"Writing out file: {out_tsv_name}", file=sys.stderr)
    with open(out_tsv_name, 'w') as out_tsv:
        # Write out a header:
        out_tsv.write(f"Gene/Transcript\tCell_Barcode\tUMI\tCount\n")

        # Write out our counts:
        for cbc in count_matrix.keys():
            for gene_tx in count_matrix[cbc]:
                for umi in count_matrix[cbc][gene_tx]:
                    out_tsv.write(
                        f"{gene_tx}\t{cbc}\t{umi}\t{count_matrix[cbc][gene_tx][umi]}\n"
                    )

    print("Done!", file=sys.stderr)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=f"Creates a count matrix TSV file from the given bam file.  "
                    f"The given bam is expected to have annotations that indicate the "
                    f"gene/transcript, Cell Barcode, and UMI for each read.  "
                    "If no Cell Barcodes are present, will quantify based on UMI only.",
        epilog="The input bam is expected to contain de-duplicated UMIs. \n"
               "The tags used by this tool are the following:\n"
               f"\tgene/transcript: {gGENE_TX_TAG_NAME}\n"
               f"\tCell Barcode:    {gCBC_TAG_NAME}\n"
               f"\tUMI:             {gUMI_TAG_NAME}\n"
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-b', '--bam',
                               help='Annotated bam file from which to create the counts matrix.',
                               required=True)
    requiredNamed.add_argument('-o', '--out-name',
                               help='Output tsv count matrix file name',
                               required=True)

    args = parser.parse_args()
    main(args.bam, args.out_name)
