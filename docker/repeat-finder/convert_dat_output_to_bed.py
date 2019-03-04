import argparse
from collections import defaultdict
from dat_utils import parse_dat_file
import tqdm

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--max-repeat-unit-size", help="Max repeat unit size (bp). Repeats with units > this length will be discarded",
        type=int, default=50)
    p.add_argument("--max-repeat-total-length", help="Max total repeat length from start to end (bp)",
        type=int, default=1000000000)
    p.add_argument("--dat-in-original-format", action="store_true", help="For .dat files generated without the trf -ngs flag")

    p.add_argument("dat_file_path")
    p.add_argument("output_file_path")
    args = p.parse_args()

    if not args.dat_file_path.endswith(".dat"):
        p.error("Invalid input filename: %s" % args.dat_file_path)

    return args


def main():

    args = parse_args()
    counters = defaultdict(int)

    with open(args.output_file_path, "w") as output_file:
        dat_records = parse_dat_file(args.dat_file_path, dat_in_original_format=args.dat_in_original_format)
        for dat_record in tqdm.tqdm(dat_records, unit=" records"):
            counters['total records'] += 1
            if "_" in dat_record.chrom or len(dat_record.chrom) >= 5:
                continue
            counters['records in chromosomes'] += 1

            if len(dat_record.repeat_unit) > args.max_repeat_unit_size:
                continue
            counters['records with repeat unit length <= %d bp' % args.max_repeat_unit_size] += 1

            if dat_record.end - dat_record.start > args.max_repeat_total_length:
                continue
            counters['records with total length <= %d bp' % args.max_repeat_total_length] += 1

            output_file.write("\t".join(map(str, [
                dat_record.chrom,
                dat_record.start,
                dat_record.end,
                dat_record.repeat_unit,
                dat_record.repeat_count,
            ])) + "\n")

    if counters['total records'] > 10 and counters['records in chromosomes'] == 0:
        raise ValueError("Chromosome parsing failed. Example chrom record: " + str(dat_record.chrom))

    for key, value in sorted(counters.items(), key=lambda x: x[1], reverse=True):
        print("%10d  %s" % (value, key))


if __name__ == "__main__":
    main()
