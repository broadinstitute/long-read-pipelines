import sys
import argparse
import numpy as np
import pysam
from functools import reduce
import operator
from tqdm import tqdm

RAW_BARCODE_TAG = 'CR'
BARCODE_QUAL_TAG = "CY"
CONF_FACTOR_SCALE = 100


def get_confidence_factor(qual_string: str, scale_factor: float = CONF_FACTOR_SCALE) -> float:
    """Get the confidence factor for the given sequence to be tallied for use with STARCODE.
    quals are assumed to be phred scale quality scores in string format and will be converted to numerical values."""
    return scale_factor * reduce(
        operator.mul, map(lambda q: 1. - 10 ** (-(ord(q) - 33.) / 10), qual_string)
    )


def main(ilmn_cellranger_bam, out_prefix):

    barcode_count_filename = f"{out_prefix}_ilmn_barcode_conf_scores.tsv"

    with open(barcode_count_filename, "a") as barcode_file:
        with pysam.AlignmentFile(ilmn_cellranger_bam, "rb", check_sq=False) as bam_file, \
                tqdm(desc=f"Processing ILMN short reads", unit="read") as pbar:
            for read in bam_file.fetch(until_eof=True):
                # This bam file should have the RAW_BARCODE_TAG and BARCODE_QUAL_TAG tags.
                # these are the data we need for our barcode file:
                barcode = read.get_tag(RAW_BARCODE_TAG)
                barcode_qual_string = read.get_tag(BARCODE_QUAL_TAG)

                # Get our conf factor and write it out to the barcode file:
                cf_raw = get_confidence_factor(barcode_qual_string)
                conf_factor = int(np.round(cf_raw))

                barcode_file.write(f"{barcode}\t{conf_factor}\n")
                pbar.update(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Process an Illumina Cellranger-processed bam file and produce a TSV file of barcodes and their "
                    "confidence scores."
    )
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-b", "--bam", help="Illumina cellranger bam file from which to extract "
                                                   "confidence scores.", required=True)
    requiredNamed.add_argument("-p", "--prefix", help='Output prefix', required=True)
    args = parser.parse_args()

    main(args.bam, args.prefix)

