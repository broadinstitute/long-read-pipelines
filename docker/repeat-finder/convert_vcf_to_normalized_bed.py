"""
This script takes a .vcf and filters it to variants where either the REF or ALT allele contains a tandem repeat.
"""

import argparse
from collections import defaultdict
import logging
import tqdm
import vcf

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger()


TRF_INPUT_SEQUENCE_MIN_LENGTH = 12


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("input_vcf_path")
    p.add_argument("output_bed_path")  # prefix for output vcf paths

    return p.parse_args()


def main(args):

    vcf_reader = vcf.VCFReader(filename=args.input_vcf_path)

    bed_file = open(args.output_bed_path, "w")

    counters = defaultdict(int)
    for i, row in tqdm.tqdm(enumerate(vcf_reader), unit=" rows"):
        if len(row.ALT) > 1:
            raise ValueError("Multi-allelic variant found. Variants should be split and normalized.")

        if row.ALT is None or row.ALT[0] is None:
            continue

        chrom, pos, ref, alt = row.CHROM, row.POS, str(row.REF), str(row.ALT[0])

        counters["total indels"] += 1
        variant_bases = ref[1:] if len(ref) > len(alt) else alt[1:]

        ins_or_del = "DEL" if len(ref) > len(alt) else "INS"
        counters[f"total {ins_or_del}"] += 1

        repeat_unit = row.INFO.get("RU", None)
        if repeat_unit is None:
            num_repeats_in_variant = 0
        else:
            num_repeats_in_variant = len(variant_bases)/float(len(repeat_unit))

            indel_size_label = "1bp" if len(repeat_unit) == 1 else (f"str" if len(repeat_unit) <= 6 else f"vntr")
            counters[f"{indel_size_label} {ins_or_del}"] += 1


        # add info about the repeat to both the ID and the INFO field for convenience & IGV
        id_field = "RU{}:{}:".format(len(repeat_unit), repeat_unit)
        if len(ref) > len(alt):
            id_field += "DEL:"
        elif len(ref) < len(alt):
            id_field += "INS:"
        id_field += "{}=>{}".format(len(ref)-1, len(alt)-1)

        bed_record = "\t".join(map(str, [chrom, pos - 1, pos - 1 + max(len(alt), len(ref)), repeat_unit, num_repeats_in_variant])) + "\n"
        bed_file.write(bed_record)
        bed_file.flush()

    bed_file.close()

    for key, value in sorted(counters.items(), key=lambda x: (x[0].split()[0], len(x[0]))):
        logger.info(f"{value:10d} {key}")


if __name__ == "__main__":
    args = parse_args()

    main(args)



