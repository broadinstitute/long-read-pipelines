#!/usr/bin/env python

import os
import sys
import argparse
from functools import reduce

import pysam
from tqdm import tqdm

gGENE_TX_TAG_NAME = "XG"
gCBC_TAG_NAME = "CB"
gUMI_TAG_NAME = "XX"


def _blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b:
            break
        yield b


def get_num_lines_in_text_file(file_name):
    with open(file_name, "r", encoding="utf-8",errors='ignore') as f:
        return sum(bl.count("\n") for bl in _blocks(f))


def main(input_bam, out_tsv_name, eq_class_tsv, gene_tx_tag, cell_barcode_tag, umi_tag):

    print("Verifying input file(s) exist...", file=sys.stderr)
    files_ok = True
    input_files = [input_bam]
    if eq_class_tsv:
        input_files.append(eq_class_tsv)
    for f in input_files:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}", file=sys.stderr)
            files_ok = False
    if not files_ok:
        sys.exit(1)
    print("Input files verified.", file=sys.stderr)

    read_to_eq_class_map = None
    if eq_class_tsv:
        print("Equivalence class TSV provided.  Ingesting read -> eq class map...", file=sys.stderr)

        read_to_eq_class_map = dict()
        print("Counting num lines in EQ class file...")
        num_lines = get_num_lines_in_text_file(f"{eq_class_tsv}")
        print(f"Num lines: {num_lines}")

        with open(eq_class_tsv, 'r') as f:
            for line in tqdm(f, desc="Reading in EQ classes", total=num_lines):
                if line.startswith("#"):
                    continue
                read_name, eq_class, gene_assignment = line.strip().split("\t")
                read_to_eq_class_map[read_name] = eq_class

    num_reads_skipped = 0
    reads_processed = 0
    reads_assigned_to_txs = 0

    print(f"Tallying gene/tx x CBC x UMI counts from {input_bam}", file=sys.stderr)
    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:

        # Create a storage system for our count matrix:
        count_matrix = dict()

        has_cbc = True

        total_reads = None
        if bam_file.has_index():
            idx_stats = bam_file.get_index_statistics()
            unaligned_reads = bam_file.nocoordinate
            aligned_reads = reduce(lambda a, b: a + b, [x.total for x in idx_stats]) if len(idx_stats) > 0 else 0
            total_reads = unaligned_reads + aligned_reads

        for read in tqdm(bam_file, desc="Progress", total=total_reads, unit=" read", file=sys.stderr):

            if read_to_eq_class_map:
                gene_tx = read_to_eq_class_map[read.query_name]
                reads_assigned_to_txs += 1
            else:
                gene_tx = read.get_tag(gene_tx_tag)

            umi = None
            try:
                umi = read.get_tag(umi_tag)
            except KeyError:
                print(f"Read {read.query_name} does not have UMI tag ({umi_tag}).", file=sys.stderr)
                num_reads_skipped += 1

            if umi:
                if has_cbc:
                    try:
                        cbc = read.get_tag(cell_barcode_tag)
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
            reads_processed += 1

    print("Done!", file=sys.stderr)
    print(f"Reads processed: {reads_processed}", file=sys.stderr)
    print(f"Reads skipped: {num_reads_skipped}", file=sys.stderr)

    print(f"Writing out file: {out_tsv_name}", file=sys.stderr)
    with open(out_tsv_name, 'w') as out_tsv:
        # Write out a header:
        if read_to_eq_class_map:
            out_tsv.write(f"Equivalence_Class\tCell_Barcode\tUMI\tCount\n")
        else:
            out_tsv.write(f"Gene/Transcript\tCell_Barcode\tUMI\tCount\n")

        # Write out our counts:
        for cbc in tqdm(count_matrix.keys(), desc="Writing gene/cbc/umi counts"):
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

    required_named_args = parser.add_argument_group('required named arguments')
    required_named_args.add_argument('-b', '--bam',
                                     help='Annotated bam file from which to create the counts matrix.',
                                     required=True)
    required_named_args.add_argument('-o', '--out-name',
                                     help='Output tsv count matrix file name',
                                     required=True)

    optional_named_args = parser.add_argument_group("optional named arguments")
    optional_named_args.add_argument('--tx-eq-class-assignments',
                                     help=f"Transcript equivalence class TSV.  If used, gene-transcript-tag will be "
                                          f"ignored and this lookup will be used.)",
                                     required=False)
    optional_named_args.add_argument('--gene-transcript-tag',
                                     type=str,
                                     help=f"Gene / Transcript BAM tag name. (Default: {gGENE_TX_TAG_NAME})",
                                     default=gGENE_TX_TAG_NAME,
                                     required=False)
    optional_named_args.add_argument('--cell-barcode-tag',
                                     type=str,
                                     help=f"Cell barcode BAM tag name. (Default: {gCBC_TAG_NAME})",
                                     default=gCBC_TAG_NAME,
                                     required=False)
    optional_named_args.add_argument('--umi-tag',
                                     type=str,
                                     help=f"UMI BAM tag name. (Default: {gUMI_TAG_NAME})",
                                     default=gUMI_TAG_NAME,
                                     required=False)

    args = parser.parse_args()
    main(args.bam, args.out_name, args.tx_eq_class_assignments, args.gene_transcript_tag, args.cell_barcode_tag,
         args.umi_tag)
