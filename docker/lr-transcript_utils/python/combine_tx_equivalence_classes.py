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

GENE_EQ_LOOKUP_FILE_SUFFIX = ".gene_equivalence_class_lookup.tsv"
GENE_EQ_DATA_FILE_SUFFIX = ".gene_name_assignments.tsv"
EQ_LOOKUP_FILE_SUFFIX = ".equivalence_class_lookup.tsv"
EQ_DATA_FILE_SUFFIX = ".equivalence_classes.tsv"

GENE_EQ_CLASS_DEFINITIONS_OUT_FILE_NAME = "gene_equivalence_class_lookup.tsv"
GENE_EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME = "gene_name_assignments.tsv"
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

gene_eq_class_def_files = dict()
gene_eq_class_data_files = dict()
tx_eq_class_def_files = dict()
tx_eq_class_data_files = dict()

if len(sys.argv) < 4:
    print("ERROR: You must provide all four equivalence class lookup files "
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
        tx_eq_class_def_files[base_name] = in_file
    elif in_file.endswith(EQ_DATA_FILE_SUFFIX):
        base_name = in_file[:in_file.find(EQ_DATA_FILE_SUFFIX)]
        tx_eq_class_data_files[base_name] = in_file
    elif in_file.endswith(GENE_EQ_LOOKUP_FILE_SUFFIX):
        base_name = in_file[:in_file.find(GENE_EQ_LOOKUP_FILE_SUFFIX)]
        gene_eq_class_data_files[base_name] = in_file
    elif in_file.endswith(GENE_EQ_DATA_FILE_SUFFIX):
        base_name = in_file[:in_file.find(GENE_EQ_DATA_FILE_SUFFIX)]
        gene_eq_class_def_files[base_name] = in_file
    else:
        print(f"ERROR: Unknown file suffix: {in_file}")
        sys.exit(1)


########################################################################################
########################################################################################
########################################################################################

def ingest_old_class_def_files(old_class_def_file_dict, file_suffix):

    old_eq_classes = dict()
    eq_class_set = set()

    for eq_def_file in sorted(list(old_class_def_file_dict.values())):
        print(f"\tProcessing {eq_def_file}")
        with open(eq_def_file, 'r') as f:
            file_key = eq_def_file[:eq_def_file.find(file_suffix)]
            old_eq_classes[file_key] = dict()
            for line in f:
                if line[0] == "#":
                    continue

                # NOTE: We don't have to reparse and sort our class definitions here because the originating
                #       script sorts them already.
                class_num, class_def = line.strip().split("\t")
                old_eq_classes[file_key][class_num] = class_def

                eq_class_set.add(class_def)

    return old_eq_classes, eq_class_set


def write_new_eq_class_def_file(out_name, eq_class_dict):
    print(f"Writing new Transcript EQ Class definitions to: {out_name}")
    with open(out_name, 'w') as f:
        f.write("#EQ_Class\tTranscript_Assignments\n")
        for eq_class, indx in eq_class_dict.items():
            f.write(f"{indx}\t{eq_class}\n")

########################################################################################
########################################################################################
########################################################################################

#########################################################
#  _____                              _       _
# |_   _| __ __ _ _ __  ___  ___ _ __(_)_ __ | |_ ___
#   | || '__/ _` | '_ \/ __|/ __| '__| | '_ \| __/ __|
#   | || | | (_| | | | \__ \ (__| |  | | |_) | |_\__ \
#   |_||_|  \__,_|_| |_|___/\___|_|  |_| .__/ \__|___/
#                                      |_|
#########################################################

# OK, we have our files collected.
# Now we need to create the new set of equivalence class labels:

print("Reading Transcript EQ Class definitions...")

old_tx_eq_classes, tx_eq_class_set = ingest_old_class_def_files(tx_eq_class_def_files, EQ_LOOKUP_FILE_SUFFIX)

# Now that we have our eq class sets, we need to renumber them:
print("Renumbering tx eq classes...")
eq_classes = {eq_class: indx for indx, eq_class in enumerate(sorted(list(tx_eq_class_set)))}

# Now we should write out our new eq class file:
write_new_eq_class_def_file(EQ_CLASS_DEFINITIONS_OUT_FILE_NAME, eq_classes)

print(f"Writing new Transcript EQ Class assignments to: {EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME}")
with open(EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME, 'w') as out_file:
    out_file.write("#Read_Name\tEQ_Class\tAssociated_Genes\n")
    for eq_class_data_file in sorted(list(tx_eq_class_data_files.values())):

        print(f"\tProcessing {eq_class_data_file}")
        file_key = eq_class_data_file[:eq_class_data_file.find(EQ_DATA_FILE_SUFFIX)]

        with open(eq_class_data_file, 'r') as f:
            for line in f:
                if line[0] == "#":
                    continue

                # Get the class info from the input file:
                fields = line.strip().split("\t")
                read_name = fields[0]
                class_num = fields[1]
                # Some eq classes don't have associated genes:
                gene_assignments = fields[2] if len(fields) > 2 else ""

                # Now lookup the new class #:
                new_class_num = eq_classes[old_tx_eq_classes[file_key][class_num]]

                # Now write out the class info:
                out_file.write(f"{read_name}\t{new_class_num}\t{gene_assignments}\n")

########################################################################################
########################################################################################
########################################################################################

############################################
#   ____
#  / ___| ___ _ __   ___  ___
# | |  _ / _ \ '_ \ / _ \/ __|
# | |_| |  __/ | | |  __/\__ \
#  \____|\___|_| |_|\___||___/
#
############################################

# OK, we have our files collected.
# Now we need to create the new set of equivalence class labels:

print("Reading Gene EQ Class definitions...")
old_gene_eq_classes, gene_eq_class_set = ingest_old_class_def_files(gene_eq_class_def_files, GENE_EQ_LOOKUP_FILE_SUFFIX)

# Now that we have our eq class sets, we need to renumber them:
print("Renumbering gene eq classes...")
gene_eq_classes = {eq_class: indx for indx, eq_class in enumerate(sorted(list(gene_eq_class_set)))}

# Now we should write out our new eq class file:
write_new_eq_class_def_file(GENE_EQ_CLASS_DEFINITIONS_OUT_FILE_NAME, gene_eq_classes)

print(f"Writing new Gene EQ Class assignments to: {GENE_EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME}")
with open(GENE_EQ_CLASS_ASSIGNMENTS_OUT_FILE_NAME, 'w') as out_file:
    for eq_class_data_file in sorted(list(gene_eq_class_data_files.values())):

        print(f"\tProcessing {eq_class_data_file}")
        file_key = eq_class_data_file[:eq_class_data_file.find(GENE_EQ_DATA_FILE_SUFFIX)]

        with open(eq_class_data_file, 'r') as f:
            for line in f:
                if line[0] == "#":
                    continue

                # Get the class info from the input file:
                fields = line.strip().split("\t")
                read_name = fields[0]
                gene_assignment = fields[1]

                # Now lookup the new class for the gene
                # Not all genes are reassigned.  In fact, most arent.
                try:
                    new_class_num = gene_eq_classes[old_gene_eq_classes[file_key][gene_assignment]]
                except KeyError:
                    new_class_num = gene_assignment

                # Now write out the class info:
                out_file.write(f"{read_name}\t{new_class_num}\n")

########################################################################################
########################################################################################
########################################################################################

end_t = time.time()
print(f"Done!  (elapsed time: {end_t - start_t}s)")
