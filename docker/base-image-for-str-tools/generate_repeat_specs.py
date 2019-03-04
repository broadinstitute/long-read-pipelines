#!/usr/bin/env python3

# ExpansionHunter requires a directory of .json files - where each .json file has the location and repeat sequence of a single STR.
# This script takes a .bed file and genome reference fasta and generates a .json file for each .bed interval
# See https://github.com/Illumina/ExpansionHunter/wiki/Inputs#repeat-specification-files  for the .json file format.

# ExpansionHunter3 requires a "variant catalog" .json file as input
# (see https://github.com/Illumina/ExpansionHunter/blob/v3.0.0-rc1/docs/04_VariantCatalogFiles.md)

# STRetch requires a .bed file where the 4th column contains the repeat unit sequence and the fourth column contains some number

import argparse
from collections import defaultdict
import gzip
import logging
import os
import pyfaidx
import unittest

logging.basicConfig(level=logging.INFO)

VALID_NUCLEOTIDES = {'A', 'C', 'G', 'T'}


def init_args():
    p = argparse.ArgumentParser()

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--expansion-hunter", action="store_true")
    g.add_argument("--expansion-hunter3", action="store_true")
    g.add_argument("--stretch", action="store_true")

    p.add_argument("-b", "--write-validated-bed", action="store_true", help="output a validated bed file")
    p.add_argument("fasta_path", help="genome reference fasta path")
    p.add_argument("bed_path", help="region path")
    p.add_argument("output_path", help="output directory path for expansion hunter, or variant_catalog.json path for expansion hunter v3")

    return p


def parse_args(argparser):
    args = argparser.parse_args()

    if not os.path.isfile(args.fasta_path):
        argparser.error("File not found: {}".format(args.fasta_path))

    if not os.path.isfile(args.bed_path):
        argparser.error("File not found: {}".format(args.bed_path))

    if args.expansion_hunter and not os.path.isdir(args.output_path):
        logging.info("Creating output dir: {}".format(args.output_path))
        os.mkdir(args.output_path)

    return args


def compute_repeat_unit(nuc_sequence):
    """Takes a nucleotide sequence and checks whether it consists of multiple exact repeats of a smaller sub-sequence.
    If yes, it returns this sub-sequence, otherwise it returns None.

    Examples:
        ACACT returns None
        ACACA returns AC
        ACTACTA returns ACT
    """
    if not nuc_sequence or len(nuc_sequence) < 2:
        raise ValueError("Invalid nuc_sequence: {}".format(nuc_sequence))

    for repeat_size in range(1, len(nuc_sequence)//2 + 1 + 1):
        repeat_unit = nuc_sequence[:repeat_size]
        repeats = repeat_unit * (1 + len(nuc_sequence)//repeat_size)
        if nuc_sequence == repeats[:len(nuc_sequence)]:
            return repeat_unit

    return None


def generate_json_file(chrom, start_1based, end_1based, repeat_unit, output_dir):
    repeat_name = "{}-{}-{}-{}".format(chrom, start_1based, end_1based, repeat_unit)

    # example: https://github.com/Illumina/ExpansionHunter/blob/master/data/repeat-specs/grch37/AR.json
    output_path = os.path.join(output_dir, "{}.json".format(repeat_name))
    #logging.info("Writing {}".format(os.path.abspath(output_path)))

    with open(output_path, "w") as f:
        f.write("""{{
    "CommonUnit": "true",
    "RepeatId": "{repeat_name}",
    "RepeatUnit": "{repeat_unit}",
    "TargetRegion": "{chrom}:{start_1based}-{end_1based}"
}}
""".format(**locals()))


def main():
    args = parse_args(init_args())
    write_validated_bed = args.write_validated_bed
    fasta_path = args.fasta_path
    bed_path = args.bed_path
    output_path = args.output_path

    # read fasta
    fasta = pyfaidx.Fasta(fasta_path, one_based_attributes=False)

    # read bed file
    fopen = gzip.open if bed_path.endswith(".gz") else open

    input_bed_file = fopen(bed_path)
    if write_validated_bed:
        output_bed_path = bed_path.replace(".bed", "_validated.bed")
        output_bed_file = fopen(output_bed_path, "w")

    counter = defaultdict(int)
    if args.expansion_hunter3 or args.stretch:
        output_file = open(output_path, "w")

    if args.expansion_hunter3:
        output_file.write("[")
        is_first_output_line = True

    for line in input_bed_file:
        if line.startswith("#"):
            continue

        counter["lines read"] += 1

        # parse .bed line
        fields = line.strip().split()

        chrom = fields[0].replace("chr", "")
        start = int(fields[1])
        end = int(fields[2])

        if "_" in chrom or len(chrom) >= 5:
            continue

        # validate bed file and repeat unit
        expected_repeat_unit = None
        reflen = None
        if len(fields) > 4 and len(set(fields[4]) - VALID_NUCLEOTIDES) == 0:
            expected_repeat_unit = fields[4]
            reflen = f"{float(end-start)/len(expected_repeat_unit):0.1f}"
        elif len(fields) > 3 and len(set(fields[3]) - VALID_NUCLEOTIDES) == 0:
            expected_repeat_unit = fields[3]
            if len(fields) > 4:
                reflen = fields[4]

        if expected_repeat_unit is not None and (end - start) % len(expected_repeat_unit) != 0:
            #logging.warn("Invalid .bed file format. 1-based start cooordinate detected: {chrom}:{start}-{end}: {expected_repeat_unit}. Coorrecting by subtracting 1 from start coord".format(**locals()))
            if not args.stretch:
                start -= 1

        nuc_sequence = str(fasta[chrom][start:end])
        repeat_unit = compute_repeat_unit(nuc_sequence)

        if repeat_unit is not None and expected_repeat_unit is not None and repeat_unit != expected_repeat_unit:
            counter["lines with repeat unit validation error"] += 1
            logging.error("Unexpected repeat sequence found at {chrom}:{start}-{end}: {repeat_unit} instead of {expected_repeat_unit}".format(**locals()))

        elif repeat_unit is None:
            counter["lines with no repeat found"] += 1
            if counter["lines with no repeat found"] < 10:
                logging.error("{chrom}:{start}-{end} - no repeat unit found in {nuc_sequence}".format(**locals()))
            elif counter["lines with no repeat found"] < 11:
                logging.error("...")

            repeat_unit = expected_repeat_unit or ""  # don't skip for stretch

        if reflen is None:
            reflen = f"{float(end-start)/len(repeat_unit)}:0.1f" if repeat_unit else 0.0

        counter["{}bp repeat unit(s)".format(len(repeat_unit))] += 1
        if args.expansion_hunter:
            generate_json_file(chrom, start+1, end, repeat_unit, output_path)
        elif args.expansion_hunter3:
            if not is_first_output_line:
                output_file.write(",")

            counter["lines written"] += 1
            output_file.write(f"""{{
    "VariantType": "Repeat",
    "LocusId": "{chrom}-{start}-{end}-{repeat_unit}",
    "LocusStructure": "({repeat_unit})*",
    "ReferenceRegion": "{chrom}:{start}-{end}"
}}""")
            is_first_output_line = False

        elif args.stretch:
            counter["lines written"] += 1
            output_file.write("\t".join(map(str, [chrom, start, end, repeat_unit, reflen])) + "\n")

        else:
            raise ValueError("Invalid state")

        if write_validated_bed:
            counter["lines written to validated bed"] += 1
            output_fields = [chrom, start, end] + ([fields[3], repeat_unit] if len(fields) > 4 else [repeat_unit])

            output_bed_file.write("\t".join(map(str, output_fields)) + "\n")

    if args.expansion_hunter3:
        output_file.write("]")

    logging.info("Parsed {} lines from {}".format(counter["lines read"], bed_path))
    for key, value in sorted(counter.items(), key=lambda x: x[1]):
        logging.info("   ==> {} {}".format(value, key))

    if args.expansion_hunter3 or args.stretch:
        logging.info(f"Done writing {counter['lines written']} lines to {output_file.name}")
        output_file.close()

    if write_validated_bed:
        output_bed_file.close()
        logging.info("Wrote {} lines to {}".format(counter["lines written to validated bed"], output_bed_path))

if __name__ == "__main__":
    main()


class Tests(unittest.TestCase):

    def test_compute_repeat_unit(self):

        self.assertEqual(compute_repeat_unit("ACA"), None)
        self.assertEqual(compute_repeat_unit("ACAG"), None)
        self.assertEqual(compute_repeat_unit("ACACA"), None)
        self.assertEqual(compute_repeat_unit("AAAA"), "A")
        self.assertEqual(compute_repeat_unit("ACAC"), "AC")
        self.assertEqual(compute_repeat_unit("CAGCAG"), "CAG")
        self.assertEqual(compute_repeat_unit("GGGCGGGC"), "GGGC")
        self.assertEqual(compute_repeat_unit("CGGGCGGG"), "CGGG")
