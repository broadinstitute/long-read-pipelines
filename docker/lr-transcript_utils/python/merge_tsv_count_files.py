#!/usr/bin/env python

# Read in data from individual quant files:
import sys
import csv
from pathlib import Path

# We must make sure we have data to merge:
if len(sys.argv) < 2:
    print(f"{sys.argv[0]} TSV_FILE_1 [TSV_FILE_2 TSV_FILE_3 ...]", file=sys.stderr)
    print(f"Merge TSV count files by naively adding their contents.", file=sys.stderr)
    print(f"Assumes all files have the same number of rows and columns and that all rows after the first row only "
          f"contain integers separated by tabs", file=sys.stderr)
    print(f"", file=sys.stderr)
    print(f"Error: you must specify at least one TSV file to merge.", file=sys.stderr)
    sys.exit(1)

file_list = sys.argv[1:]
MERGED_OUT_TSV = "merged_counts.tsv"

# Make sure all files exist:
for f in file_list:
    p = Path(f)
    if not p.is_file():
        print(f"Error: TSV file does not exist: {f}", file=sys.stderr)
        sys.exit(1)

# Get the header:
print(f"Getting header...", file=sys.stderr)
with open(file_list[0], 'r') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter="\t")
    for row in tsvreader:
        header_fields = row
        break

# Merge the contents:
print(f"Merging contents...", file=sys.stderr)
out_rows = []
is_first_file = True
for f in file_list:
    with open(f, 'r') as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t")
        # Skip the first row:
        next(tsvreader)
        for i, row in enumerate(tsvreader):
            int_row = [int(c) for c in row]
            if is_first_file:
                out_rows.append(int_row)
            else:
                tmp_row = int_row
                for j in range(len(int_row)):
                    tmp_row[j] += out_rows[i][j]
                out_rows[i] = tmp_row

        is_first_file = False

print(f"Writing output...", file=sys.stderr)
with open(MERGED_OUT_TSV, 'w') as out_file:
    # Write out header:
    out_file.write("\t".join(header_fields))
    out_file.write("\n")

    for row in out_rows:
        out_row = '\t'.join([str(x) for x in row])
        out_file.write(f"{out_row}\n")

print(f"Done!", file=sys.stderr)
