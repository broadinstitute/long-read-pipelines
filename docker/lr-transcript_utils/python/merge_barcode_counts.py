#!/usr/bin/env python

# Read in data from individual quant files:
import sys
import csv
from pathlib import Path
from tqdm import tqdm

# We must make sure we have data to merge:
if len(sys.argv) != 2:
    print(f"{sys.argv[0]} BARCODE_TSV_FILE", file=sys.stderr)
    print(f"Merge barcode counts in a given TSV count file.", file=sys.stderr)
    print(f"TSV file must be unheadered and must have 2 columns: BARCODE COUNT.", file=sys.stderr)
    print(f"All barcodes with the same name will have their counts added.", file=sys.stderr)
    print(f"Resulting TSV file will have unique barcodes and consiladated counts.", file=sys.stderr)
    sys.exit(1)

tsv_file_name = sys.argv[1]
MERGED_OUT_TSV = "merged_counts.tsv"

# Make sure all files exist:
p = Path(tsv_file_name)
if not p.is_file():
    print(f"Error: TSV file does not exist: {tsv_file_name}", file=sys.stderr)
    sys.exit(1)

with open(tsv_file_name, 'r') as f, tqdm(desc="Processing Barcodes", unit=" barcode") as pbar:
    barcode_count_dict = dict()
    tsvreader = csv.reader(f, delimiter="\t")
    for barcode, count in tsvreader:
        count = int(count)
        try:
            barcode_count_dict[barcode] += count
        except KeyError:
            barcode_count_dict[barcode] = count
        pbar.update(1)

print(f"Num Barcodes: {len(barcode_count_dict)}", file=sys.stderr)

print(f"Writing output to {MERGED_OUT_TSV}", file=sys.stderr)
with open(MERGED_OUT_TSV, "w") as f, \
        tqdm(desc="Writing barcode file", unit=" barcode", total=len(barcode_count_dict)) as pbar:
    for barcode, count in barcode_count_dict.items():
        f.write(f"{barcode}\t{count}\n")
        pbar.update(1)

print(f"Done!", file=sys.stderr)
