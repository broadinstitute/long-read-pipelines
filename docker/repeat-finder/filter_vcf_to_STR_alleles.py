"""
This script takes a .vcf and filters it to variants where either the REF or ALT allele contains a tandem repeat.
"""

import argparse
from collections import defaultdict
import logging
import numpy as np
import os
import pyfaidx
from scipy.stats import entropy
import tempfile
import tqdm
import unittest
import vcf
from dat_utils import TRFRunner

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger()


TRF_INPUT_SEQUENCE_MIN_LENGTH = 12


def parse_args():
    p = argparse.ArgumentParser()

    # VCF processing params
    p.add_argument("-R", "--reference-fasta-path", help="Reference genome fasta path", required=True)

    # TRF params
    p.add_argument("--trf-input-sequence-min-length", type=int, default=TRF_INPUT_SEQUENCE_MIN_LENGTH)
    p.add_argument("--trf-path", help="Tandem Repeats Finder command path", default="/trf409.linux")
    p.add_argument("--mismatch-penalty", type=int, default=7)
    p.add_argument("--indel-penalty", type=int, default=7)
    p.add_argument("--minscore", type=int, default=24)

    # post-filter
    p.add_argument("--min-flanking-repeats-in-reference", type=int, default=4, help="If a repeat is found in the variant, "
        "require the flanking reference sequence to have at least this many additional copies of the repeat for the  variant "
        "to be considered a TR. Otherwise, the variant will be discarded.")
    p.add_argument("--min-flanking-sequence-size", type=int, default=0, help="How many base pairs of flanking sequence to "
        "consider in the above filter. If not specified or 0, the value will be the length of the repeat unit in the variant * min_flanking_repeats_in_reference + 1")

    p.add_argument("input_vcf_path")
    p.add_argument("output_prefix")  # prefix for output vcf paths

    return p.parse_args()


def compute_entropy(nuc_sequence):
    value, counts = np.unique(list(nuc_sequence), return_counts=True)
    return entropy(counts)


def get_flanking_reference_sequences(fasta_obj, chrom, pos, ref, alt, num_flanking_bases=None):
    """Takes a single insertion or deletion and returns flanking reference sequence around the variant.

    Args:
        fasta_obj: pyfaidx wrapper of the reference genome fasta
        chrom (str): variant chromosome
        pos (str): variant position
        ref (str): variant ref allele
        alt (str): variant alt allele
        num_flanking_bases (int): (optional) num flanking bases to return. If not specified, it will be set to 5x the
            allele length.

    Return:
        3-tuple: left_flanking_reference_sequence,
            variant_bases (the bases that were inserted or deleted by this variant),
            right_flanking_reference_sequence
    """
    if len(ref) == 1 and len(alt) > 1:
        variant_bases = alt[1:]
    elif len(ref) > 1 and len(alt) == 1:
        variant_bases = ref[1:]
    else:
        raise ValueError(f"Invalid variant: {chrom}-{pos}-{ref}-{alt}. Make sure all multi-allelic variants have been split and normalized.")

    if num_flanking_bases is None:
        num_flanking_bases = 5 * len(variant_bases)

    left_start = max(0, pos - num_flanking_bases)
    right_end = min(pos + num_flanking_bases, len(fasta_obj[chrom]))

    left_flanking_reference_sequence = fasta_obj[chrom][left_start : pos]
    right_flanking_reference_sequence = fasta_obj[chrom][pos : right_end]

    return left_flanking_reference_sequence, variant_bases, right_flanking_reference_sequence


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


INDEL_SIZE_LABELS = ["all_sizes"]   # ["1bp", "str", "vntr", "all_sizes"]

def main(args):

    vcf_reader = vcf.VCFReader(filename=args.input_vcf_path)

    vcf_writers = {}
    #bed_files = {}
    for indel_size_label in INDEL_SIZE_LABELS:
         for ins_or_del in ["DEL", "INS", "INDEL"]:
             #bed_files[(ins_or_del, indel_size_label)] = open(f"{args.output_prefix}_{ins_or_del}_{indel_size_label}.bed", "w")
             vcf_writers[(ins_or_del, indel_size_label)] = vcf.VCFWriter(open(f"{args.output_prefix}_{ins_or_del}_{indel_size_label}.vcf", "w"), vcf_reader)

    vcf_denovo_writer = vcf.VCFWriter(open(f"{args.output_prefix}_denovo.vcf", "w"), vcf_reader)

    trf_runner = TRFRunner(args.trf_path, mismatch_penalty=args.mismatch_penalty, indel_penalty=args.indel_penalty, minscore=args.minscore)

    def reverse_string(s):
        """Utility function"""
        if s is None:
            return s

        assert isinstance(s, str)

        return s[::-1]

    def find_repeat_unit(nucleotide_sequence, fraction_covered_by_repeat=None, only_return_this_sorted_repeat_unit=None, verbose=False):
        """Takes a nucleotide sequence and looks for tandem repeats - using simple exact matching if nucleotide_sequence is very short,
        or running Tandem Repeat Finder if it's longer.

        Args:
            nucleotide_sequence (str): short tandem repeat sequence
            fraction_covered_by_repeat (float): (optional) what fraction of the nucleotide sequence must be covered by a repeat motif
                before the repeat is returned
            only_return_this_sorted_repeat_unit (str): if specified, a repeat unit will only be returned if it has this sequence when sorted in alphabetical order
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

            if only_return_this_sorted_repeat_unit:
                dat_records = [dr for dr in dat_records if "".join(sorted(dr.repeat_unit)) == only_return_this_sorted_repeat_unit]

            if not dat_records:
                return None, 0

            dat_records = sorted(dat_records, key=lambda dat_record: (dat_record.alignment_score, dat_record.repeat_count), reverse=True)
            repeat_unit, repeat_count = dat_records[0].repeat_unit, dat_records[0].repeat_count

        if repeat_unit is None:
            return None, 0

        if only_return_this_sorted_repeat_unit is not None and "".join(sorted(repeat_unit)) != only_return_this_sorted_repeat_unit:
            return None, 0

        return repeat_unit, repeat_count

    counters = defaultdict(int)
    fasta_obj = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False, as_raw=True)
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
            counters[f"skipping len(repeat_unit) == 1 and num_repeats_in_variant == 1"] += 1
            continue

        counters[f"total {ins_or_del} with repeat"] += 1
        indel_size_label = "1bp" if len(repeat_unit) == 1 else (f"str" if len(repeat_unit) <= 6 else f"vntr")
        counters[f"{indel_size_label} {ins_or_del} with repeat"] += 1

        left_flanking_reference_sequence, variant_bases, right_flanking_reference_sequence = get_flanking_reference_sequences(
            fasta_obj, chrom, pos, ref, alt, num_flanking_bases=max(args.min_flanking_sequence_size, len(repeat_unit) * (args.min_flanking_repeats_in_reference + 1)))

        #right_repeat_unit, num_repeats_right = find_repeat_unit(right_flanking_reference_sequence)
        #left_repeat_unit, num_repeats_left = find_repeat_unit(left_flanking_reference_sequence)
        #if (left_repeat_unit != repeat_unit or num_repeats_left < 4) and \
        #   (right_repeat_unit != repeat_unit or num_repeats_right < 4):
        #    continue

        sorted_repeat_unit = "".join(sorted(repeat_unit))

        # for the left side, reverse and then reverse back so that find_repeat_unit starts from the variant bases
        # rather than from the end of the flanking sequence. That way, the length of the flanking sequence doesn't matter as much.
        left_sequence = reverse_string(left_flanking_reference_sequence + variant_bases)
        left_repeat_unit, num_repeats_left = find_repeat_unit(left_sequence, only_return_this_sorted_repeat_unit=sorted_repeat_unit)
        left_repeat_unit = reverse_string(left_repeat_unit)

        right_sequence = variant_bases + right_flanking_reference_sequence
        right_repeat_unit, num_repeats_right = find_repeat_unit(right_sequence, only_return_this_sorted_repeat_unit=sorted_repeat_unit)

        sorted_right_repeat_unit = "".join(sorted(right_repeat_unit)) if right_repeat_unit is not None else None
        right_match = sorted_right_repeat_unit == sorted_repeat_unit and num_repeats_right >= num_repeats_in_variant + args.min_flanking_repeats_in_reference
        left_match = False

        left_entropy = compute_entropy(left_flanking_reference_sequence)/len(left_flanking_reference_sequence)
        right_entropy = compute_entropy(right_flanking_reference_sequence)/len(right_flanking_reference_sequence)

        is_denovo = False
        if not right_match:
            sorted_left_repeat_unit = "".join(sorted(left_repeat_unit)) if left_repeat_unit is not None else None
            left_match = sorted_left_repeat_unit == sorted_repeat_unit and num_repeats_left >= num_repeats_in_variant + args.min_flanking_repeats_in_reference
            if not left_match:
                if num_repeats_in_variant > 4 and len(repeat_unit) > 1:
                    new_INFO = {
                        'RU': repeat_unit,
                        "LeftEntropy": str(left_entropy),
                        "RightEntropy": str(right_entropy),
                    }
                    print("\ndenovo repeat unit: ", repeat_unit)
                    print("variant: ", variant_bases)
                    print("\nleft repeat unit: ", left_repeat_unit, "  entropy: ", left_entropy, "  left_sequence:", left_flanking_reference_sequence)
                    print("\nright repeat unit: ", right_repeat_unit, "  entropy: ", right_entropy, "  right_sequence:", right_flanking_reference_sequence)

                    is_denovo = True
                    record = vcf.model._Record(
                        chrom,
                        pos,
                        repeat_unit,
                        row.REF,
                        row.ALT,
                        row.QUAL,
                        row.FILTER,
                        new_INFO,
                        row.FORMAT,
                        sample_indexes={c.sample: c for c in row.samples},
                        samples=row.samples)

                    vcf_denovo_writer.write_record(record)

                else:
                    if len(repeat_unit) > 1:
                        counters[f"skipping {ins_or_del} with only {num_repeats_in_variant} repeat in variant"] += 1

                    continue

        #bed_start = ((pos - len(left_flanking_reference_sequence)) if left_match else pos)
        #bed_end = (pos + len(right_flanking_reference_sequence)) if right_match else pos
        #while tuple(sorted(str(fasta_obj[chrom][bed_end : bed_end + len(repeat_unit)]))) == sorted_repeat_unit:
        #    bed_end += len(repeat_unit)
        #    continue

        num_repeats = num_repeats_left if left_match else num_repeats_right

        #print(f"------ Not Denovo: {chrom} {pos} {repeat_unit} {num_repeats_in_variant}", chrom, bed_start, bed_end)

        counters[f"total alleles with repeat and {'right' if right_repeat_unit == repeat_unit else 'left'} flanking repeat"] += 1
        counters[f"{indel_size_label} {ins_or_del} with repeat and {'right' if right_repeat_unit == repeat_unit else 'left'} flanking repeat"] += 1

        # add info about the repeat to both the ID and the INFO field for convenience & IGV
        id_field = "RU{}:{}:".format(len(repeat_unit), repeat_unit)
        if len(ref) > len(alt):
            id_field += "DEL:"
        elif len(ref) < len(alt):
            id_field += "INS:"
        id_field += "{}=>{}".format(len(ref)-1, len(alt)-1)

        new_INFO = dict(row.INFO)
        new_INFO.update({
            'RU': repeat_unit,
            "RepeatClass": indel_size_label,
            "RepeatCountInVariant": str(num_repeats_in_variant),
            "RepeatCountInVariantAndRef": str(max(num_repeats_left, num_repeats_right, num_repeats_left + num_repeats_right - num_repeats_in_variant)),
            "RepeatCountInVariantAndRefLeft": str(num_repeats_left),
            "RepeatCountInVariantAndRefRight": str(num_repeats_right),
            "LacksMatchingRepeatsInRef": str(is_denovo),
            "LeftEntropy": str(left_entropy),
            "RightEntropy": str(right_entropy),
            "LeftFlank": str(left_flanking_reference_sequence),
            "RightFlank": str(right_flanking_reference_sequence),
        })

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

        #bed_record = "\t".join(map(str, [chrom, bed_start, bed_end, repeat_unit, num_repeats])) + "\n"

        for ins_or_del_i in [ins_or_del, "INDEL"]:
            for indel_size_label_j in INDEL_SIZE_LABELS:  # [indel_size_label, "all_sizes"]
                vcf_writers[(ins_or_del_i, indel_size_label_j)].write_record(record)
                vcf_writers[(ins_or_del_i, indel_size_label_j)].flush()
                #bed_files[(ins_or_del_i, indel_size_label_j)].write(bed_record)
                #bed_files[(ins_or_del_i, indel_size_label_j)].flush()

    for ins_or_del in ["DEL", "INS", "INDEL"]:
        for indel_size_label in INDEL_SIZE_LABELS:
            vcf_writers[(ins_or_del, indel_size_label)].close()
            #bed_files[(ins_or_del, indel_size_label)].close()

    for key, value in sorted(counters.items(), key=lambda x: (x[0].split()[0], len(x[0]))):
        logger.info(f"{value:10d} {key}")


if __name__ == "__main__":
    args = parse_args()

    main(args)




class Tests(unittest.TestCase):
    """

        Insertion:

             1234567890--1234567890..
        REF: TGTGTGTGTT--ACACACACTGTGTGTGTGGGGGG

        ALT: TGTGTGTGTTacACACACACTGTGTGTGTGGGGGG

            output-left: TGTGTGTTac
            output-right: acACACACAC

            pos: 10
            ref: T
            alt: TAC


        Deletion:

             12345678901234567890
        REF: TGTGTGTGTTACACACACTGTGTGTGTGGGGGG

        ALT: TGTGTGTGTT--ACACACTGTGTGTGTGGGG

            output-left: TGTGTGTTAC
            output-right: ACACACACTG

            pos: 10
            ref: TAC
            alt: T
    """

    def setUp(self):
        self.temp_fasta_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
        self.temp_fasta_file.write(">chrTest1\n")
        self.temp_fasta_file.write("TGTGTGTGTTACACACACTGTGTGTGTGGGGGG\n")


        self.temp_fasta_file.write(">chrTest2\n")
        self.temp_fasta_file.write("TGTGTGTGTTACACACTGTGTGTGTGGGGGG\n")


        self.test3_seq1 = "GTG"*5
        self.test3_seq2 = "ATC"*5
        self.test3_seq3 = "TTC"*5
        self.temp_fasta_file.write(">chrTest3\n")
        self.temp_fasta_file.write(f"{self.test3_seq1}{self.test3_seq2}{self.test3_seq3}\n")
        self.temp_fasta_file.close()

        self.fasta_obj = pyfaidx.Fasta(self.temp_fasta_file.name, one_based_attributes=False, as_raw=True)

    def test_get_containing_sequences1(self):
        # insertion 2-STR
        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest1", 10, "T", "TAC")
        self.assertEqual(left_sequence, "TGTGTGTGTT")
        self.assertEqual(variant_sequence, "AC")
        self.assertEqual(right_sequence, "ACACACACTG")

        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest1", 10, "T", "TACAC")
        self.assertEqual(left_sequence, "TGTGTGTGTT")
        self.assertEqual(variant_sequence, "ACAC")
        self.assertEqual(right_sequence, "ACACACACTGTGTGTGTGGG")

        # deletion
        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest2", 10, "TAC", "T")
        self.assertEqual(left_sequence, "TGTGTGTGTT")
        self.assertEqual(variant_sequence, "AC")
        self.assertEqual(right_sequence, "ACACACTGTG")

        # insertion 3-STR
        left_sequence, variant_sequence, right_sequence = get_flanking_reference_sequences(self.fasta_obj, "chrTest3", 15, self.test3_seq1[-1], self.test3_seq1[-1] + self.test3_seq2)
        self.assertEqual(left_sequence, self.test3_seq1)
        self.assertEqual(variant_sequence, self.test3_seq2)
        self.assertEqual(right_sequence, self.test3_seq2 + self.test3_seq3)


    def test_find_exact_repeat_unit(self):
        self.assertEqual(find_exact_repeat_unit("C"), (None, 0))
        self.assertEqual(find_exact_repeat_unit("AA"), ("A", 2))
        self.assertEqual(find_exact_repeat_unit("AC"), (None, 0))
        self.assertEqual(find_exact_repeat_unit("ACC"), (None, 0))
        self.assertEqual(find_exact_repeat_unit("ACAG"), (None, 0))
        self.assertEqual(find_exact_repeat_unit("ACTAG"), (None, 0))
        self.assertEqual(find_exact_repeat_unit("ACTATAGATGATAG"), (None, 0))
        self.assertEqual(find_exact_repeat_unit("CCCCG"), (None, 0))

        self.assertEqual(find_exact_repeat_unit("ACA"), (None, 0))
        self.assertEqual(find_exact_repeat_unit("ACACA"), ("AC", 2))
        self.assertEqual(find_exact_repeat_unit("ACTACTA"), ("ACT", 2))
        self.assertEqual(find_exact_repeat_unit("AAAA"), ("A", 4))
        self.assertEqual(find_exact_repeat_unit("ACAC"), ("AC", 2))
        self.assertEqual(find_exact_repeat_unit("CAGCAG"), ("CAG", 2))
        self.assertEqual(find_exact_repeat_unit("GGGCGGGC"), ("GGGC", 2))
        self.assertEqual(find_exact_repeat_unit("CGGGCGGG"), ("CGGG", 2))


        trf_runner = TRFRunner(os.path.expanduser("~/bin/trf409.macosx"), mismatch_penalty=3, indel_penalty=5, minscore=8)

        repeat_sequence = "ACACACACAC"
        dat_records = list(trf_runner.run_TRF(repeat_sequence))
        self.assertEqual(len(dat_records), 1)
        self.assertEqual(dat_records[0].repeat_unit, "AC")
        self.assertEqual(dat_records[0].repeat_count, 5)
        self.assertEqual(find_exact_repeat_unit(repeat_sequence), ("AC", 5))

        repeat_sequence = "ACGAACGAACGAACGAATGA"
        dat_records = list(trf_runner.run_TRF(repeat_sequence))
        self.assertEqual(len(dat_records), 1)
        self.assertEqual(dat_records[0].repeat_unit, "ACGA")
        self.assertEqual(dat_records[0].repeat_count, 5)
        self.assertEqual(find_exact_repeat_unit(repeat_sequence), (None, 0))

    def tearDown(self):
        if os.path.isfile(self.temp_fasta_file.name):
            os.remove(self.temp_fasta_file.name)
