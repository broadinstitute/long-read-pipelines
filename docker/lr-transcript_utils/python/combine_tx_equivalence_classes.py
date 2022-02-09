#!/usr/bin/env python

# This is a quick script to combine all the equivalence class data from multiple disjoint subsets of reads into a final
# combined set of files - one eq class lookup and one set of eq classes.
#
# This is based on output assumptions from the `quantify_gff_reads.py` script (listed explicitly below).
#
# Author: Jonn Smith
# Date:  2022 02 09

import sys
import time

EQ_LOOKUP_FILE_SUFFIX = ".equivalence_class_lookup.tsv"
EQ_DATA_FILE_SUFFIX = ".equivalence_classes.tsv"

EQ_CLASS_DEFINITIONS_OUT_FILE_NAME = "equivalence_class_lookup.tsv"
EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME = "equivalence_classes.tsv"

# No real arg parsing here.
# We take a collection of files and split them into two types:
#     1) eq definitions (lookup tables)
#     2) data files themselves.

# Assumptions:
# 1) Each EQ Data file has a corresponding Lookup file
# 2) The names of the EQ data files and Lookup files have the same prefixes up until the suffixes defined above.
# 3) The contents of each EQ class in the definitions file (i.e. the transcripts and class codes) are sorted.

start_t = time.time()

eq_class_def_files = dict()
eq_class_data_files = dict()

if len(sys.argv) < 2:
    print("ERROR: You must provide both equivalence class lookup files "
          "and equivalence class assignment files.", file=sys.stderr)
    sys.exit(1)

if (len(sys.argv) - 1) % 2 != 0:
    print("ERROR: An uneven number of files was provided!  "
          "At least one file doesn't have corresponding lookups/assignments.", file=sys.stderr)
    sys.exit(1)

print("Classifying input files...")
for in_file in sys.argv[1:]:
    if in_file.endswith(EQ_LOOKUP_FILE_SUFFIX):
        base_name = in_file[:in_file.find(EQ_LOOKUP_FILE_SUFFIX)]
        eq_class_def_files[base_name] = in_file
    elif in_file.endswith(EQ_DATA_FILE_SUFFIX):
        base_name = in_file[:in_file.find(EQ_DATA_FILE_SUFFIX)]
        eq_class_data_files[base_name] = in_file
    else:
        print(f"ERROR: Unknown file suffix: {in_file}")
        sys.exit(1)

# OK, we have our files collected.
# Now we need to create the new set of equivalence class labels:
old_eq_classes = dict()
eq_class_set = set()

print("Reading EQ Class definitions...")
for eq_def_file in sorted(list(eq_class_def_files.values())):
    print(f"\tProcessing {eq_def_file}")
    with open(eq_def_file, 'r') as f:
        file_key = eq_def_file[:eq_def_file.find(EQ_LOOKUP_FILE_SUFFIX)]
        old_eq_classes[file_key] = dict()
        for line in f:
            if line[0] == "#":
                continue

            # NOTE: We don't have to reparse and sort our class definitions here because the originating
            #       script sorts them already.
            class_num, class_def = line.strip().split("\t")
            class_num = int(class_num)
            old_eq_classes[file_key][class_num] = class_def

            eq_class_set.add(class_def)

# Now that we have our eq class sets, we need to renumber them:
print("Renumbering eq classes...")
eq_classes = {eq_class: indx for indx, eq_class in enumerate(sorted(list(eq_class_set)))}

# Now we should write out our new eq class file:
print(f"Writing new EQ Class definitions to: {EQ_CLASS_DEFINITIONS_OUT_FILE_NAME}")
with open(EQ_CLASS_DEFINITIONS_OUT_FILE_NAME, 'w') as f:
    f.write("#EQ_Class\tTranscript_Assignments\n")
    for eq_class, indx in eq_classes.items():
        f.write(f"{indx}\t{eq_class}\n")

print(f"Writing new EQ Class assignments to: {EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME}")
with open(EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME, 'w') as out_file:
    out_file.write("#Read_Name\tEQ_Class\tAssociated_Genes\n")
    for eq_class_data_file in sorted(list(eq_class_data_files.values())):

        print(f"\tProcessing {eq_class_data_file}")
        file_key = eq_class_data_file[:eq_class_data_file.find(EQ_DATA_FILE_SUFFIX)]

        with open(eq_class_data_file, 'r') as f:
            for line in f:
                if line[0] == "#":
                    continue

                # Get the class info from the input file:
                fields = line.strip().split("\t")
                read_name = fields[0]
                class_num = int(fields[1])
                # Some eq classes don't have associated genes:
                gene_assignments = fields[2] if len(fields) > 2 else ""

                # Now lookup the new class #:
                new_class_num = eq_classes[old_eq_classes[file_key][class_num]]

                # Now write out the class info:
                out_file.write(f"{read_name}\t{new_class_num}\t{gene_assignments}\n")

end_t = time.time()
print(f"Done!  (elapsed time: {end_t - start_t}s)")
