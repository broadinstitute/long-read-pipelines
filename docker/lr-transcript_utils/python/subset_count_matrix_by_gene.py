#!/usr/bin/env python

# Read in data from individual quant files:
import sys
import glob
import csv
import time
import os

import anndata
import scanpy as sc

# Right now we accept only 1 argument - the directory in which the quant files live:
if len(sys.argv) < 3:
    print(f"{sys.argv[0]} COUNTS_MATRIX_TSV GENENAME1 [GENENAME2 GENENAME3 ...]", file=sys.stderr)
    print(f"Error: you must specify at least one gene name.", file=sys.stderr) 
    sys.exit(1)

count_matrix = sys.argv[1]
gene_names = set([f"|{g}|" for g in sys.argv[2:]])

base_name = "count_matrix_subset_by_gene"
tsv_name = f"{base_name}.tsv"

print(f"Count matrix file: {count_matrix}", file=sys.stderr)
print(f"Selected genes: {', '.join(sys.argv[2:])}", file=sys.stderr)

print(f"Getting gene-specific header...", file=sys.stderr)
with open(count_matrix, 'r') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter="\t") 
    for row in tsvreader:
        field_names = row
        break
fields_to_keep = set()
for i,field in enumerate(field_names):
    for g in gene_names:
        if g in field:
            fields_to_keep.add(i)
            break
new_header = [field_names[i] for i in fields_to_keep]
print(f"Keeping fields: ", file=sys.stderr)
for i,n in zip(fields_to_keep, new_header):
    print(f"    {i}: {n}", file=sys.stderr)

if len(fields_to_keep) == 0:
    print("Error: Did not keep any data.  Did you name your genes correctly (check gencode)?", file=sys.stderr)
    sys.exit(1)

print(f"Creating TSV: {tsv_name}", file=sys.stderr)
with open(tsv_name, 'w') as out_file:
    # Write out header:
    out_file.write("\t")
    out_file.write("\t".join(new_header))
    out_file.write("\n")

    with open(count_matrix, 'r') as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t") 
        # Skip the first row:
        next(tsvreader)
        for row in tsvreader:
            out_row = '\t'.join([row[i] for i in fields_to_keep])
            out_file.write(f"{row[0]}\t{out_row}\n")

print('Loading in TSV file...', file=sys.stderr)
adata = sc.read(tsv_name, ext='tsv')
print('Saving anndata...', file=sys.stderr)
adata.write(f"{base_name}.h5ad")
print('Done!', file=sys.stderr)

