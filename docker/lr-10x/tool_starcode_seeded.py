import subprocess

import pysam
import argparse
from Bio.Seq import Seq
import time
import os.path
import gzip
import ctypes
import sys
import shutil

import numpy as np

from tqdm import tqdm
from functools import reduce
import operator
import array

sys.path.append('/lrma')
from ssw import ssw_lib

BARCODE_TAG = 'CB'
RAW_BARCODE_TAG = 'CR'

BARCODE_CORRECTED_TAG = "XC"


def read_barcodes(barcodes_filename):
    """
    Reads a line-separated file of barcodes. If a line ends with '-1', this suffix is removed. A gzip file is supported if the extension of the file is .gz
    :param barcodes_filename: Filename of the barcode file.
    :return: A set of barcodes.
    """
    barcodes = set()
    _, barcodes_file_extension = os.path.splitext(barcodes_filename)
    is_gzip = barcodes_file_extension == '.gz'
    with (open(barcodes_filename) if not is_gzip else gzip.open(barcodes_filename)) as barcodes_file:
        for line in barcodes_file:
            if is_gzip:
                line = line.decode('utf-8')
            if line[-3:-1] == '-1':
                barcodes.add(line[:-3])
            else:
                barcodes.add(line[:-1])
    return barcodes


def _process_starcode_stdout(starcode_proc, starcode_output_file = None):
    num_clusters = 0
    correction_dict = dict()

    for cluster_line in starcode_proc.stdout:
        num_clusters += 1
        cluster_line = cluster_line.decode("ascii")
        cluster_data = cluster_line.strip().split('\t')
        cluster = cluster_data[0]
        cluster_sequences = cluster_data[2].split(',')
        for cluster_sequence in cluster_sequences:
            correction_dict[cluster_sequence] = cluster

        if starcode_output_file is not None:
            starcode_output_file.write(cluster_line)

    return correction_dict, num_clusters


def perform_barcode_correction_starcode_with_barcode_counts(
        barcode_count_filename, analysis_name, starcode_path,
        starcode_dist, starcode_cluster_ratio, starcode_sphere, starcode_connected_comp
):
    """
    Performs the barcode correction using starcode
    :param analysis_name: Prefix for the stats files
    :param starcode_path: Relative or absolute path to the starcode executable
    :param barcode_count_filename: The tsv file containing barcode counts.  If not None, this file will be used by
    starcode to create the correction dictionary.
    :return: A dictionary with the raw barcodes as keys and the corresponding corrected barcodes as values
    """

    starcode_args = [starcode_path, "--print-clusters", "--quiet", "-i", barcode_count_filename]
    if starcode_dist:
        starcode_args.append("--dist")
        starcode_args.append(starcode_dist)
    if starcode_cluster_ratio:
        starcode_args.append("--cluster-ratio")
        starcode_args.append(starcode_cluster_ratio)
    if starcode_sphere:
        starcode_args.append("--sphere")
        starcode_args.append(starcode_sphere)
    if starcode_connected_comp:
        starcode_args.append("--connected-comp")
        starcode_args.append(starcode_connected_comp)

    print(f"Executing starcode as: {' '.join(starcode_args)}", file=sys.stderr)
    starcode_proc = subprocess.Popen(starcode_args, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    starcode_proc.stdin.close()

    starcode_output_file = open(analysis_name + '_starcode.tsv', 'w') if analysis_name is not None else None
    try:
        if starcode_output_file is not None:
            starcode_output_file.write('cluster\toccurrence\tindices\n')

        correction_dict, num_clusters = _process_starcode_stdout(starcode_proc, starcode_output_file)

    finally:
        if starcode_output_file is not None:
            starcode_output_file.close()

    return correction_dict


def main(bam_filename, counts_filename, analysis_name, whitelist_10x_filename, starcode_path,
         starcode_dist, starcode_cluster_ratio, starcode_sphere, starcode_connected_comp):
    """
    Main function for the tool
    """

    if whitelist_10x_filename:
        print(f"Reading in 10x whitelist: {whitelist_10x_filename}")
        whitelist_10x = read_barcodes(whitelist_10x_filename)
    else:
        whitelist_10x = None

    print('Performing barcode corrections...')
    correction_dict = perform_barcode_correction_starcode_with_barcode_counts(
        counts_filename, analysis_name, starcode_path,
        starcode_dist, starcode_cluster_ratio, starcode_sphere, starcode_connected_comp
    )

    correction_dict['.'] = '.'
    print(f"Barcode Correction Dict Length: {len(correction_dict)}")

    with pysam.AlignmentFile(bam_filename, 'rb', check_sq=False) as bam_file, \
            tqdm(desc=f"Correcting barcodes", unit="read") as pbar:

        with pysam.AlignmentFile(analysis_name + '_corrected_barcodes.bam', 'wb', check_sq=False, header=bam_file.header) as output_file:
            for read in bam_file.fetch(until_eof=True):

                # Get our raw barcode and set some basic outputs:
                raw_barcode = read.get_tag(RAW_BARCODE_TAG)
                read.set_tag(BARCODE_TAG, raw_barcode, value_type='Z')
                read.set_tag(BARCODE_CORRECTED_TAG, False)

                # Get the corrected output:
                corrected_barcode = correction_dict[raw_barcode]

                # If we used a 10x whitelist, only set the corrected barcode if it's in the whitelist:
                if whitelist_10x is not None:
                    if corrected_barcode in whitelist_10x:
                        read.set_tag(BARCODE_TAG, corrected_barcode, value_type='Z')
                        read.set_tag(BARCODE_CORRECTED_TAG, True)
                else:
                    read.set_tag(BARCODE_TAG, corrected_barcode, value_type='Z')
                    read.set_tag(BARCODE_CORRECTED_TAG, True)

                output_file.write(read)
                pbar.update(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads an input BAM file which has been annotated with a raw cell barcode and corrects the existing'
                    'cell barcode using a counts file and starcode.'
                    'The results are written out to another BAM file.')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-b', '--bam', help='BAM filename', required=True)
    requiredNamed.add_argument('-c', '--counts', help='Counts TSV filename', required=True)
    requiredNamed.add_argument('-n', '--name', help='Analysis name (output prefix)', required=True)

    parser.add_argument('--whitelist-10x', help='10x whitelist filename. This may be GZIP compressed (has to have extension .gz in that case)')
    parser.add_argument('--whitelist-illumina', help='Illumina whitelist filename. This may be GZIP compressed (has to have extension .gz in that case)')
    parser.add_argument('--starcode-path', help='Path to the starcode executable', type=str, default='/lrma/starcode-master/starcode')

    # Options for starcode:
    optionalStarcode = parser.add_argument_group("Starcode Options")
    optionalStarcode.add_argument("--dist", help="Defines the maximum Levenshtein distance for clustering.", required=False)
    optionalStarcode.add_argument("--cluster-ratio", help="(Message passing only) Specifies the minimum sequence count ratio to cluster two matching sequences.",
                                  required=False)
    optionalStarcode.add_argument("--sphere", help="Use sphere clustering algorithm instead of message passing (MP).", action="store_true")
    optionalStarcode.add_argument("--connected-comp", help="Clusters are defined by the connected components.", action="store_true")

  # general options:
  #   -d --dist: maximum Levenshtein distance (default auto)
  #   -t --threads: number of concurrent threads (default 1)
  #   -q --quiet: quiet output (default verbose)
  #   -v --version: display version and exit
  #
  # cluster options: (default algorithm: message passing)
  #   -r --cluster-ratio: min size ratio for merging clusters in
  #              message passing (default 5.0)
  #   -s --sphere: use sphere clustering algorithm
  #   -c --connected-comp: cluster connected components

    args = parser.parse_args()

    # Process args for mutex options:
    # Some starcode options are mutually exclusive:
    if args.cluster_ratio and (args.sphere or args.connected_comp):
        raise RuntimeError("Arguments are mutually exclusive: `--cluster-ratio`, `--sphere`, `--connected-comp`.")

    if args.spheres and (args.cluster_ratio or args.connected_comp):
        raise RuntimeError("Arguments are mutually exclusive: `--cluster-ratio`, `--sphere`, `--connected-comp`.")

    if args.connected_comp and (args.sphere or args.cluster_ratio):
        raise RuntimeError("Arguments are mutually exclusive: `--cluster-ratio`, `--sphere`, `--connected-comp`.")

    if args.whitelist_illumina and not args.whitelist_10x:
        print('Illumina whitelist provided but no 10x whitelist provided.')
        exit(1)

    main(args.bam, args.counts, args.name, args.whitelist_10x, args.starcode_path, args.dist, args.cluster_ratio, args.sphere, args.connected_comp)

