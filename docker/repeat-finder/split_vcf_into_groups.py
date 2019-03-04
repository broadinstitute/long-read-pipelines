"""
This script takes a .vcf and filters it to variants where either the REF or ALT allele contains a tandem repeat.
"""

import argparse
from collections import defaultdict
import logging
import tqdm
import vcf
from dat_utils import TRFRunner

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger()


TRF_INPUT_SEQUENCE_MIN_LENGTH = 12


def parse_args():
    p = argparse.ArgumentParser()

    # TRF params
    p.add_argument("--trf-input-sequence-min-length", type=int, default=TRF_INPUT_SEQUENCE_MIN_LENGTH)
    p.add_argument("--trf-path", help="Tandem Repeats Finder command path", default="/trf409.linux")
    p.add_argument("--mismatch-penalty", type=int, default=7)
    p.add_argument("--indel-penalty", type=int, default=7)
    p.add_argument("--minscore", type=int, default=24)

    p.add_argument("input_vcf_path")
    p.add_argument("output_prefix")  # prefix for output vcf paths

    return p.parse_args()


def find_exact_repeat_unit(nuc_sequence):
    """Takes a nucleotide sequence and checks if it consists of multiple exact repeats.
    If yes, it returns the repeat unit and number of repeats, otherwise it returns (None, 0).

    Examples:
        ACACT returns (None, 0)
        ACACA returns (AC, 2)
        ACTACTA returns (ACT, 2)
    """
    if len(nuc_sequence) > 1:
        for repeat_size in range(1, (len(nuc_sequence) + 1)//2 + 1):
            repeat_unit = nuc_sequence[:repeat_size]
            repeats = repeat_unit * (1 + len(nuc_sequence)//repeat_size)
            repeats_count = len(nuc_sequence)//len(repeat_unit)
            if nuc_sequence == repeats[:len(nuc_sequence)] and repeats_count >= 2:
                return repeat_unit, repeats_count

    return None, 0


def main(args):

    vcf_reader = vcf.VCFReader(filename=args.input_vcf_path)

    vcf_writers = {}
    for indel_size_label in ["1bp", "str", "vntr", "all_sizes"]:
         for ins_or_del in ["DEL", "INS", "INDEL"]:
             vcf_writers[(ins_or_del, indel_size_label)] = vcf.VCFWriter(open(f"{args.output_prefix}_{ins_or_del}_{indel_size_label}.vcf", "w"), vcf_reader)

    trf_runner = TRFRunner(args.trf_path, mismatch_penalty=args.mismatch_penalty, indel_penalty=args.indel_penalty, minscore=args.minscore)

    def find_repeat_unit(nucleotide_sequence, fraction_covered_by_repeat=None, only_return_this_repeat_unit=None, verbose=False):
        """Takes a nucleotide sequence and looks for tandem repeats - using simple exact matching if nucleotide_sequence is very short,
        or running Tandem Repeat Finder if it's longer.

        Args:
            nucleotide_sequence (str): short tandem repeat sequence
            fraction_covered_by_repeat (float): (optional) what fraction of the nucleotide sequence must be covered by a repeat motif
                before the repeat is returned
            only_return_this_repeat_unit (str): if specified, a repeat unit will only be returned if it has this sequence
            verbose (bool): Print debug messages.

        Returns:
            2-tuple: repeat-unit string, number of repeats found in nucleotide_sequence
        """
        if len(nucleotide_sequence) == 1:
            repeat_unit, repeat_count = nucleotide_sequence, 1

        elif len(nucleotide_sequence) < args.trf_input_sequence_min_length:
            repeat_unit, repeat_count = find_exact_repeat_unit(nucleotide_sequence)

        else:
            dat_records = list(trf_runner.run_TRF(nucleotide_sequence))

            if fraction_covered_by_repeat is not None:
                dat_records = [dr for dr in dat_records
                    if dr.repeat_count * len(dr.repeat_unit) / float(len(nucleotide_sequence)) > fraction_covered_by_repeat]

            if only_return_this_repeat_unit:
                dat_records = [dr for dr in dat_records if dr.repeat_unit == only_return_this_repeat_unit]

            if not dat_records:
                return None, 0

            dat_records = sorted(dat_records, key=lambda dat_record: (dat_record.alignment_score, dat_record.repeat_count), reverse=True)
            repeat_unit, repeat_count = dat_records[0].repeat_unit, dat_records[0].repeat_count

        if only_return_this_repeat_unit is not None and repeat_unit != only_return_this_repeat_unit:
            return None, 0

        return repeat_unit, repeat_count

    counters = defaultdict(int)
    for i, row in tqdm.tqdm(enumerate(vcf_reader), unit=" rows"):
        if len(row.ALT) > 1:
            raise ValueError("Multi-allelic variant found. Variants should be split and normalized.")

        if row.ALT is None or row.ALT[0] is None:
            continue

        chrom, pos, ref, alt = row.CHROM, row.POS, str(row.REF), str(row.ALT[0])

        counters["all variants"] += 1
        if len(ref) == len(alt):
            continue

        counters["total indels"] += 1
        variant_bases = ref[1:] if len(ref) > len(alt) else alt[1:]

        ins_or_del = "DEL" if len(ref) > len(alt) else "INS"
        counters[f"total {ins_or_del}"] += 1

        repeat_unit, num_repeats_in_variant = find_repeat_unit(variant_bases, fraction_covered_by_repeat=0.90)
        if repeat_unit is None:
            repeat_unit = variant_bases
            num_repeats_in_variant = 1
            #continue

        if len(repeat_unit) == 1 and num_repeats_in_variant == 1:
            continue

        counters[f"total {ins_or_del} with repeat"] += 1
        indel_size_label = "1bp" if len(repeat_unit) == 1 else (f"str" if len(repeat_unit) <= 6 else f"vntr")
        counters[f"{indel_size_label} {ins_or_del} with repeat"] += 1

        # add info about the repeat to both the ID and the INFO field for convenience & IGV
        id_field = "RU{}:{}:".format(len(repeat_unit), repeat_unit)
        if len(ref) > len(alt):
            id_field += "DEL:"
        elif len(ref) < len(alt):
            id_field += "INS:"
        id_field += "{}=>{}".format(len(ref)-1, len(alt)-1)

        new_INFO = dict(row.INFO)
        new_INFO['RU'] = repeat_unit

        record = vcf.model._Record(
            chrom,
            pos,
            id_field,
            row.REF,
            row.ALT,
            row.QUAL,
            row.FILTER,
            new_INFO,
            row.FORMAT,
            sample_indexes={c.sample: c for c in row.samples},
            samples=row.samples)

        for ins_or_del_i in [ins_or_del, "INDEL"]:
            for indel_size_label_j in [indel_size_label, "all_sizes"]:
                vcf_writers[(ins_or_del_i, indel_size_label_j)].write_record(record)

    for indel_size_label in ["1bp", "str", "vntr", "all_sizes"]:
        for ins_or_del in ["DEL", "INS", "INDEL"]:
            vcf_writers[(ins_or_del, indel_size_label)].close()

    for key, value in sorted(counters.items(), key=lambda x: (x[0].split()[0], len(x[0]))):
        logger.info(f"{value:10d} {key}")


if __name__ == "__main__":
    args = parse_args()

    main(args)

