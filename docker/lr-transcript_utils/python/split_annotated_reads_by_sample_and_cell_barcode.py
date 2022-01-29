#!/usr/bin/env python

import argparse
import pysam


UNKNOWN_READ_GROUP_NAME = 'UNKNOWN'


def write_manifest(manifest_cb_fasta_dict_dict, outname='flair_reads_manifest.tsv'):
    """
    Write out a manifest file to disk for use by FLAIR.

    A manifest file is of the form:
    sample1	conditionA	batch1	./sample1_reads.fq
    sample1	conditionA	batch2	./sample1_2_reads.fq
    sample2	conditionA	batch1	./sample2_reads.fq
    ...
    sampleN	conditionZ	batchM	./sampleN_reads.fq

    :param manifest_cb_fasta_dict_dict: Dict of sample_name -> {CELL_BARCODE, FASTA_FILE}
    :param outname: File name of the resulting manifest file.
    """

    default_condition = 'A'

    with open(outname, 'w') as f:
        for sample, barcode_fasta_dict in manifest_cb_fasta_dict_dict.items():
            for barcode, fasta in barcode_fasta_dict.items():
                f.write(f'{sample}\t{default_condition}\t{barcode}\t{fasta}\n')


def main(bam_filename, out_base_name):
    """
    Do the work!

    :param bam_filename: Filename of the reads BAM file
    :param out_base_name: Output base filename.
    """

    # Dict of sample_name -> {CELL_BARCODE, FASTA_FILE}
    manifest_cb_fasta_dict_dict = dict()

    tot_count = 0
    split_count = 0
    update_interval = 10000

    # Now go through the BAM file and add the read group annotations:
    with pysam.AlignmentFile(bam_filename, 'rb', check_sq=False) as bam_file:

        for read in bam_file.fetch(until_eof=True):

            split_count += 1
            tot_count += 1

            sample_name = read.get_tag('RG')

            # Skip unknowns since we can't pin down their origin:
            if sample_name == UNKNOWN_READ_GROUP_NAME:
                continue

            # Get the cell barcode:
            if read.has_tag('CB'):
                cell_barcode = read.get_tag('CB')
            elif read.has_tag('CR'):
                cell_barcode = read.get_tag('CR')
            else:
                # no CR / CB annotation means we can't determine the barcode and we must skip:
                continue

            fasta_name = f'{out_base_name}_{sample_name}_{cell_barcode}.fasta'
            cb_fasta_dict = {cell_barcode: fasta_name}

            if sample_name in manifest_cb_fasta_dict_dict:
                if cell_barcode in manifest_cb_fasta_dict_dict[sample_name]:
                    with open(fasta_name, 'a') as f:
                        f.write(f'>{read.query_name}\n')
                        f.write(f'{read.query_sequence}\n')
                else:
                    with open(fasta_name, 'w') as f:
                        f.write(f'>{read.query_name}\n')
                        f.write(f'{read.query_sequence}\n')

                    manifest_cb_fasta_dict_dict[sample_name][cell_barcode] = fasta_name
            else:
                with open(fasta_name, 'w') as f:
                    f.write(f'>{read.query_name}\n')
                    f.write(f'{read.query_sequence}\n')

                manifest_cb_fasta_dict_dict[sample_name] = cb_fasta_dict

            if split_count % update_interval == 0:
                print(f"Reads split:\t{split_count}\t[{tot_count}]")

    write_manifest(manifest_cb_fasta_dict_dict, out_base_name + '_flair_reads_manifest.tsv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Splits an annotated bam into individual FASTA files, each containing reads from a '
                    'single sample + cell barcode pair.'
    )
    requiredNamed = parser.add_argument_group('Required Named Arguments')
    requiredNamed.add_argument('-b', '--bam', help='BAM filename', required=True)
    requiredNamed.add_argument('-o', '--outbasename', help='Output base filename', required=True)

    args = parser.parse_args()

    main(args.bam, args.outbasename)

