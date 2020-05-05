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
if len(sys.argv) < 2:
    print(f"Error: you must specify a directory containing quant files produced by salmon.", file=sys.stderr) 
    sys.exit(1)

d = sys.argv[1]

# Some basic argument validation:
if not os.path.isdir(d):
    print(f"Error: Given directory does not exist: {d}", file=sys.stderr)
    sys.exit(1)

print(f"Getting list of quant.sf files from {d} ...", file=sys.stderr)
file_list = glob.glob(os.path.join(d, '*quant.sf'))

base_name = "differential_expression-known_isoforms"
tsv_name = f"{base_name}.tsv"

print(f"Detected Files: {', '.join(file_list)}")

print("Getting list transcripts ...", file=sys.stderr)
var_names = []
field_names = []
with open(file_list[0], 'r') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter="\t") 
    is_first = True
    for row in tsvreader:
        if is_first:
            field_names = row
            is_first = False
            continue
        var_names.append(row[0])

print(f"Creating TSV: {tsv_name}", file=sys.stderr)
with open(tsv_name, 'w') as out_file:
    # Write out header:
    out_file.write("\t")
    out_file.write("\t".join(var_names))
    out_file.write("\n")

    print(f"Reading in sample data ({len(file_list)} files)...", file=sys.stderr)
    start_t = time.time()
    cur_t = start_t
    i = 0
    step_height = 50
    for f in file_list:
        sample_name = f[-30:-9]
        out_file.write(f"{sample_name}")
        with open(f, 'r') as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t") 
            next(tsvreader)
            for row in tsvreader:
                # We actually only need the NumReads data (5th column)
                out_file.write(f"\t{row[4]}")
        out_file.write("\n")
    
        i += 1
        if i % step_height == 0:
            l_cur_t = cur_t
            cur_t = time.time()
            elapsed_t = cur_t - start_t
            remaining_t = (i / elapsed_t / step_height) * (len(file_list) - i)
            print(f"    Files Read:\t{i}/{len(file_list)} ({i/len(file_list)*100.1:.2f}%) - dt: {cur_t - l_cur_t:.2f} - Elapsed Time: {elapsed_t:.2f} - Remaining Time: {remaining_t:.2f}s", file=sys.stderr)

#        if i == 100:
#            print("DEBUG: Terminating early.", file=sys.stderr)
#            break

print('Loading in TSV file...', file=sys.stderr)
adata = sc.read(tsv_name, ext='tsv').transpose()
print('Saving anndata...', file=sys.stderr)
adata.write(f"{base_name}.h5ad")
print('Done!', file=sys.stderr)
