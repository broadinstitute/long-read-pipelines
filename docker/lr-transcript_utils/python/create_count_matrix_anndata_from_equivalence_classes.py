#!/usr/bin/env python
import math
import os
import sys
import argparse
import csv
import pickle

import scipy
import anndata
import numpy as np
import pandas as pd

from tqdm import tqdm

GENE_ENTRY_STRING = "gene"
TX_ENTRY_STRING = "transcript"

CONTIG_FIELD = "contig"
START_FIELD = "start"
END_FIELD = "end"

GENE_ID_FIELD = "gene_id"
TX_ID_FIELD = "transcript_id"

GENCODE_GENE_NAME_FIELD = "gene_name"
GENCODE_TX_NAME_FIELD = "transcript_name"

STRINGTIE_GENE_ID_FIELD = "ref_gene_id"
STRINGTIE_TX_ID_FIELD = "reference_id"
STRINGTIE_GENE_NAME_FIELD = "ref_gene_name"

ALT_NAME_SUFFIX = "_PAR_Y"


def get_gtf_field_val_dict(gtf_file, entry_type_filter=TX_ENTRY_STRING, force_rebuild=False):

    global TX_ENTRY_STRING
    global GENCODE_GENE_NAME_FIELD
    global GENE_ID_FIELD
    global GENCODE_TX_NAME_FIELD
    global TX_ID_FIELD
    global ALT_NAME_SUFFIX

    global CONTIG_FIELD
    global START_FIELD
    global END_FIELD

    gtf_dict = dict()

    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file does not exist: {gtf_file}")

    pickle_file_name = os.path.splitext(os.path.basename(gtf_file))[0] + f".{entry_type_filter}.big_map.pickle"

    if not force_rebuild and os.path.exists(pickle_file_name):
        print(f"Loading GTF field map from {pickle_file_name}...", end="\t", file=sys.stderr)
        gtf_dict = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        print(f"Creating master GTF field map keyed by transcript ID using {gtf_file}...",
              flush=True, file=sys.stderr)

        with open(gtf_file, "r") as f, tqdm(desc="Processing GTF File", unit=" line") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            for row in tsv_file:
                # Ignore comments and make sure we only process transcript entries:
                if not row[0].startswith("#") and row[2] == entry_type_filter:

                    # Parse the row data into a dict and make sure we're only looking at protein coding rows:
                    row_data_dict = {
                        field.strip().split(" ")[0].replace('"', ""): field.strip().split(" ")[1].replace('"', "") for
                        field in row[8].split(";") if len(field) != 0
                    }
                    row_data_dict[CONTIG_FIELD] = row[0]
                    row_data_dict[START_FIELD] = int(row[3])
                    row_data_dict[END_FIELD] = int(row[4])

                    # Make sure our names are unique:
                    try:
                        if row_data_dict[TX_ID_FIELD].endswith(ALT_NAME_SUFFIX):
                            row_data_dict[GENCODE_TX_NAME_FIELD] = row_data_dict[GENCODE_TX_NAME_FIELD] + ALT_NAME_SUFFIX
                            row_data_dict[GENCODE_GENE_NAME_FIELD] = row_data_dict[GENCODE_GENE_NAME_FIELD] + ALT_NAME_SUFFIX
                    except KeyError:
                        if row_data_dict[GENE_ID_FIELD].endswith(ALT_NAME_SUFFIX):
                            row_data_dict[GENCODE_GENE_NAME_FIELD] = row_data_dict[GENCODE_GENE_NAME_FIELD] + ALT_NAME_SUFFIX

                    # Add this row to our dict keyed by transcript ID:
                    try:
                        gtf_dict[row_data_dict[TX_ID_FIELD]] = row_data_dict
                    except KeyError:
                        gtf_dict[row_data_dict[GENE_ID_FIELD]] = row_data_dict

                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(gtf_dict, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    return gtf_dict


def intervals_overlap(contig, start, end, contig2, start2, end2):
    if contig == contig2:
        if (start <= start2 <= end) or (start <= end2 <= end) or (start2 <= start and end2 >= end):
            return True
    return False


def interval_overlaps_any_in_interval_list(contig, start, end, interval_list):
    for contig2, start2, end2 in interval_list:
        if intervals_overlap(contig, start, end, contig2, start2, end2):
            return True
    return False


def get_approximate_gencode_gene_assignments(gtf_field_dict, gencode_field_val_dict, overlap_threshold=0.1):
    gene_name_assignments = []
    gene_id_assignments = []
    ambiguity_markers = np.empty(len(gtf_field_dict), dtype=bool)

    # Make structures here that will allow fast search:
    gencode_contigs = np.array([f[CONTIG_FIELD] for f in gencode_field_val_dict.values()])
    gencode_starts = np.array([f[START_FIELD] for f in gencode_field_val_dict.values()])
    gencode_ends = np.array([f[END_FIELD] for f in gencode_field_val_dict.values()])

    gencode_index_name_map = {i: k for i, k in enumerate(gencode_field_val_dict.keys())}

    with tqdm(desc="Processing GTF Entries", unit=" entry", total=len(gtf_field_dict)) as pbar:
        for i, (k, v) in enumerate(gtf_field_dict.items()):
            contig = v[CONTIG_FIELD]
            start = v[START_FIELD]
            end = v[END_FIELD]

            # now we test for inclusion:
            gencode_overlapping_indices = np.where(
                (contig == gencode_contigs) &
                (((gencode_starts <= start) & (start <= gencode_ends)) |
                 ((gencode_starts <= end) & (end <= gencode_ends)) |
                 ((start <= gencode_starts) & (gencode_ends <= end)))
            )[0]

            # If we have some overlaps, then we need to mark them:
            if np.any(gencode_overlapping_indices):

                max_gencode_index = 0
                overlap_scores = np.zeros(len(gencode_overlapping_indices))
                for j, overlap_index in enumerate(gencode_overlapping_indices):
                    # Determine the amount of overlap in the two ranges:
                    overlap_start = max(start, gencode_starts[overlap_index])
                    overlap_end = min(end, gencode_ends[overlap_index])
                    overlap_scores[j] = (overlap_end - overlap_start)

                    # Store the max here to make it a little faster:
                    if overlap_scores[j] > overlap_scores[max_gencode_index]:
                        max_gencode_index = j

                is_ambiguous = (min(overlap_scores) / max(overlap_scores) > overlap_threshold) if len(gencode_overlapping_indices) > 1 else False

                ####################################
                # DEBUGGING:
                print(f"Genes overlapping [{i}:{k} @ {contig}:{start}-{end}] [ambiguous: {is_ambiguous}]:")
                for j in range(len(gencode_overlapping_indices)):
                    key = gencode_index_name_map[gencode_overlapping_indices[j]]
                    c = gencode_field_val_dict[key][CONTIG_FIELD]
                    s = gencode_field_val_dict[key][START_FIELD]
                    e = gencode_field_val_dict[key][END_FIELD]

                    best_string = "***" * 2 if j == max_gencode_index else ""

                    print(f"\t{j}\t{gencode_overlapping_indices[j]}\t{key}\t{gencode_field_val_dict[key][GENCODE_GENE_NAME_FIELD]} @ {c}:{s}-{e} ({overlap_scores[j]}){best_string}")
                ####################################

                # Set our gene as the one with the most overlap:
                key = gencode_index_name_map[gencode_overlapping_indices[max_gencode_index]]
                gene_name_assignments.append(gencode_field_val_dict[key][GENCODE_GENE_NAME_FIELD])
                gene_id_assignments.append(gencode_field_val_dict[key][GENE_ID_FIELD])
                ambiguity_markers[i] = is_ambiguous
            else:
                # We have no existing transcripts for which to add annotations.
                # We must add the label of the de-novo gene name and mark as unambiguous:
                gene_name_assignments.append(v[GENE_ID_FIELD])
                gene_id_assignments.append(v[GENE_ID_FIELD])
                ambiguity_markers[i] = False

            print(f"Assigned: {k} -> {gene_name_assignments[i]}<{gene_id_assignments[i]}>  (ambiguous: {ambiguity_markers[i]})")
            pbar.update(1)

    gene_name_assignments = np.array(gene_name_assignments)
    gene_id_assignments = np.array(gene_id_assignments)

    return gene_name_assignments, gene_id_assignments, ambiguity_markers


def create_combined_anndata(input_tsv, tx_eq_class_def_map, gene_eq_class_def_map, read_eq_class_map,
                            gtf_field_dict, gencode_gtf_field_dict, overlapping_gene_name_set=None,
                            overlap_intervals_label="overlaps_intervals_of_interest",
                            force_recount=False):

    """Create an anndata object holding the given gene/transcript information.
    NOTE: This MUST be a sparse matrix - we got lots of data here.

    First convert the counts to a matrix that looks like what scanpy expects:
    columns = Transcripts name / Transcript ID / Gene Name / Gene ID (variables)
    Rows = cell barcodes (observations)
    data = counts"""

    global GENCODE_GENE_NAME_FIELD
    global GENE_ID_FIELD
    global GENCODE_TX_NAME_FIELD
    global TX_ID_FIELD

    global STRINGTIE_GENE_ID_FIELD
    global STRINGTIE_TX_ID_FIELD
    global STRINGTIE_GENE_NAME_FIELD

    global CONTIG_FIELD
    global START_FIELD
    global END_FIELD

    cell_barcode_to_tx_eq_class_to_umi_dict = dict()
    cell_barcode_to_tx_eq_class_count_dict = dict()

    pickle_file_name = os.path.splitext(os.path.basename(input_tsv))[0] + ".tx_raw_count_matrix.pickle"

    if not force_recount and os.path.exists(pickle_file_name):
        print(f"Loading cell count map from {pickle_file_name}...", end="\t", file=sys.stderr)
        cell_barcode_to_tx_eq_class_count_dict = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        # Get our cell tx counts:
        with open(input_tsv, "r") as f, tqdm(desc="Processing Raw Cell Counts", unit=" count") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            next(tsv_file)
            for row in tsv_file:
                tx_eq_class_id = row[0]
                cell_barcode = row[1]
                umi = row[2]

                # Handle the Transcript names:
                if cell_barcode not in cell_barcode_to_tx_eq_class_to_umi_dict:
                    cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode] = dict()
                    cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode][tx_eq_class_id] = set()
                    cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode][tx_eq_class_id].add(umi)

                    cell_barcode_to_tx_eq_class_count_dict[cell_barcode] = dict()
                    cell_barcode_to_tx_eq_class_count_dict[cell_barcode][tx_eq_class_id] = 1

                elif tx_eq_class_id not in cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode]:
                    cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode][tx_eq_class_id] = set()
                    cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode][tx_eq_class_id].add(umi)

                    cell_barcode_to_tx_eq_class_count_dict[cell_barcode][tx_eq_class_id] = 1

                elif umi not in cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode][tx_eq_class_id]:
                    cell_barcode_to_tx_eq_class_to_umi_dict[cell_barcode][tx_eq_class_id].add(umi)

                    cell_barcode_to_tx_eq_class_count_dict[cell_barcode][tx_eq_class_id] += 1

                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(cell_barcode_to_tx_eq_class_count_dict, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    ##################################################################################################################
    # Write our cell TX counts to our adata object:

    # Create unique row / column identifiers into which to aggregate data:
    cell_barcodes = np.array(list(cell_barcode_to_tx_eq_class_count_dict.keys()))
    tx_eq_classes = np.unique(np.array(list(tx_eq_class_def_map.keys())))

    tx_eq_class_index_dict = {name: i for i, name in enumerate(tx_eq_classes)}

    # Populate the count matrix:
    pickle_file_name = os.path.splitext(os.path.basename(input_tsv))[0] + ".cell_transcript_count_matrix.pickle"
    if not force_recount and os.path.exists(pickle_file_name):
        print(f"Loading count map from {pickle_file_name}...", end="\t", file=sys.stderr)
        count_mat = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        count_mat = scipy.sparse.lil_matrix((len(cell_barcodes), len(tx_eq_classes)), dtype=np.uint32)
        with tqdm(desc=f"Creating cell transcript eq class count matrix", unit=" cell",
                  total=len(cell_barcode_to_tx_eq_class_count_dict)) as pbar:

            for i, (cb, counts_dict) in enumerate(cell_barcode_to_tx_eq_class_count_dict.items()):
                # Put the counts for each transcript in the right indices:
                for tx_eq_class, count in counts_dict.items():
                    count_mat[i, tx_eq_class_index_dict[tx_eq_class]] = count
                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(count_mat, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    ############################################################################################################
    # Now we set up the variables that we're going to apply to each observation:
    # We'll try two different ways - one for Gencode GTFs and one for stringtie GTFs:

    # First we handle our overlapping tests:
    # Mark the equivalence classes that overlap our intervals:
    tx_eq_class_overlap_flags = None
    if overlapping_gene_name_set:
        print(f"Determining overlap intervals for the given interval list... ", end="", file=sys.stderr)

        tx_eq_class_overlap_flags = np.empty(len(tx_eq_classes), dtype=bool)

        seen_eq_classes = set()
        for read_name, (tx_eq_class, gene_eq_class) in read_eq_class_map.items():

            # We only need to do this once per each eq class:
            if tx_eq_class in seen_eq_classes:
                continue

            try:
                gene_ids = gene_eq_class_def_map[gene_eq_class]
            except KeyError:
                # The gene is not in the eq class map and so we have just this one
                gene_ids = tuple([gene_eq_class])

            does_overlap = False
            for gene in gene_ids:
                # Remove the version number from the gene if it's present:
                idx = gene.find(".")
                if idx != -1:
                    gene = gene[:gene.find(".")]

                does_overlap = gene in overlapping_gene_name_set
                if does_overlap:
                    break

            if does_overlap:
                tx_eq_class_overlap_flags[tx_eq_class_index_dict[tx_eq_class]] = True

            # Make sure we don't see this guy again:
            seen_eq_classes.add(tx_eq_class)

        print("Done!", file=sys.stderr)

    # Let's create a lookup from geneID -> gene name:
    gene_name_lookup_by_id = dict()
    for _, gencode_gene_info in gencode_gtf_field_dict.items():
        gene_name_lookup_by_id[gencode_gene_info[GENE_ID_FIELD]] = gencode_gene_info[GENCODE_GENE_NAME_FIELD]

    # Now we handle the rest of our metadata (var / obs info):
    gene_eq_classes = [""] * len(tx_eq_classes)
    transcript_ids = [""] * len(tx_eq_classes)
    gene_ids = [""] * len(tx_eq_classes)
    gene_names = [""] * len(tx_eq_classes)

    is_de_novo = np.empty(len(tx_eq_classes), dtype=bool)
    is_gene_id_ambiguous = np.empty(len(tx_eq_classes), dtype=bool)

    seen_eq_classes = set()
    for read_name, (tx_eq_class, gene_eq_class) in read_eq_class_map.items():

        if tx_eq_class in seen_eq_classes:
            continue

        tx_indx = tx_eq_class_index_dict[tx_eq_class]

        # Get our gene ids:
        try:
            read_gene_ids = sorted(gene_eq_class_def_map[gene_eq_class])
        except KeyError:
            # The gene is not in the eq class map and so we have just this one
            read_gene_ids = tuple([gene_eq_class])

        tx_ids = sorted(tx_eq_class_def_map[tx_eq_class])

        # Set gene eq class:
        gene_eq_classes[tx_indx] = gene_eq_class

        # Set transcript ids:
        transcript_ids[tx_indx] = ",".join([txid for txid, cc in tx_ids])

        # Set gene IDs:
        gene_ids[tx_indx] = ",".join(read_gene_ids)

        # Get the gene names:
        if gene_eq_class.startswith("ENSG"):
            gene_name_assignment = gene_name_lookup_by_id[gene_eq_class]
        else:
            assigned_gene_names = []
            for gid in read_gene_ids:
                try:
                    assigned_gene_names.append(gene_name_lookup_by_id[gid])
                except KeyError:
                    assigned_gene_names.append(gid)

            gene_name_assignment = ",".join(assigned_gene_names)

        gene_names[tx_indx] = gene_name_assignment

        # Any transcript with an equivalence class containing a stringtie 2 transcript, or the null transcript is de novo:
        is_de_novo[tx_indx] = any([txid.startswith("STRG") or txid.startswith("-") for txid, cc in tx_ids])

        # Any gene is ambiguous if it has a gene equivalence class:
        # NOTE: This might have to change by the multi-mapping gene / read pairs.
        is_gene_id_ambiguous = len(read_gene_ids) != 1

        seen_eq_classes.add(tx_eq_class)

    # Create our anndata object now:
    count_adata = anndata.AnnData(count_mat.tocsr())

    # Add our variables:
    col_df = pd.DataFrame()
    col_df["transcript_eq_classes"] = tx_eq_classes
    col_df["gene_eq_classes"] = gene_eq_classes
    col_df["transcript_ids"] = transcript_ids
    col_df["gene_ids"] = gene_ids

    col_df["gene_names"] = gene_names

    col_df["is_de_novo"] = is_de_novo
    col_df["is_gene_id_ambiguous"] = is_gene_id_ambiguous

    # If we're doing interval overlaps, add our label:
    if overlapping_gene_name_set:
        print(f"Adding {overlap_intervals_label} column for overlapping intervals...", file=sys.stderr, end="")
        col_df[f"{overlap_intervals_label}"] = tx_eq_class_overlap_flags
        print("Done!", file=sys.stderr)

    # Assign the data to our anndata object:
    count_adata.var = col_df
    count_adata.var_names = [str(t) for t in tx_eq_classes]

    # Add our observations:
    row_df = pd.DataFrame()
    row_df["Cell Barcode"] = cell_barcodes

    count_adata.obs = row_df
    count_adata.obs_names = cell_barcodes

    return count_adata


def read_intervals_from_tsv(filename):
    with open(filename, "r") as f, tqdm(desc="Processing intervals file", unit=" line") as pbar:
        intervals = []
        tsv_file = csv.reader(f, delimiter="\t")
        for row in tsv_file:
            if (not row[0].startswith("#")) and (not row[0].startswith("@")):
                try:
                    contig = row[0]
                    start = int(row[1])
                    end = int(row[2])
                    intervals.append((contig, start, end))
                except IndexError:
                    print("ERROR: Input interval file does not seem to be a properly formatted TSV.  "
                          "It must have the format: CONTIG    START    END", file=sys.stderr)
                    sys.exit(1)
            pbar.update(1)
    return intervals


def read_gene_names_from_intervals_file(filename):
    # Intervals file is of the format:
    # File provenance: https://www.genenames.org/data/genegroup/#!/group/370 , Gencode 37
    # # HGNC:12281	TRGJP2	T cell receptor gamma joining P2	Approved	T cell receptor gene	TCRGJP2	JP2	7p14.1	6972	ENSG00000211688	OTTHUMG00000155222	375	T cell receptor gamma locus at 7p14
    # chr14 21797287  21797886
    # ...
    #
    # So we can parse the gene info lines and get the gene names for overlaps instead of the intervals!
    with open(filename, "r") as f, tqdm(desc="Processing intervals file", unit=" line") as pbar:
        gene_names = set()
        tsv_file = csv.reader(f, delimiter="\t")
        for row in tsv_file:
            if row[0].startswith("# File provenance"):
                continue
            if not row[0].startswith("#"):
                continue
            gene_name = row[10]
            gene_names.add(gene_name)
            pbar.update(1)
    return gene_names


def parse_tx_eq_class_defs_file(file_name):
    eq_class_map = dict()
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            eq_class, raw_class_assignments = line.strip().split("\t")
            class_assignments = []
            for raw_class in raw_class_assignments.split(","):
                eq_class_name, cc = raw_class.split(";")
                class_assignments.append((eq_class_name, cc))
            eq_class_map[eq_class] = tuple(class_assignments)

    return eq_class_map


def parse_gene_eq_class_assignments_file(file_name):
    eq_class_map = dict()
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            eq_class_name, gene_assignments = line.strip().split("\t")
            gene_assignments = gene_assignments.split(',')
            eq_class_map[eq_class_name] = tuple(gene_assignments)

    return eq_class_map


def parse_tx_eq_class_assignments(file_name):
    read_eq_class_map = dict()
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            read_name, tx_eq_class, gene_assignment = line.strip().split("\t")
            read_eq_class_map[read_name] = tuple([tx_eq_class, gene_assignment])

    return read_eq_class_map


def main(input_tsv, gtf_file, out_prefix,
         tx_eq_class_definitions,
         tx_eq_class_assignments,
         gene_eq_class_definitions,
         gene_eq_class_assignments,
         overlap_interval_filename=None, overlap_intervals_label=None,
         gencode_reference_gtf=None):

    print("Verifying input file(s) exist...", file=sys.stderr)
    files_ok = True
    files_to_check = [input_tsv, gtf_file, tx_eq_class_definitions, tx_eq_class_assignments,
                      gene_eq_class_definitions, gene_eq_class_assignments]
    if overlap_interval_filename:
        files_to_check.append(overlap_interval_filename)
    if gencode_reference_gtf:
        files_to_check.append(gencode_reference_gtf)
    for f in files_to_check:
        if not os.path.exists(f):
            print(f"ERROR: Input file does not exist: {f}", file=sys.stderr)
            files_ok = False
    if not files_ok:
        sys.exit(1)

    overlapping_gene_names = None
    if overlap_interval_filename:
        print("Verifying contents of {overlap_interval_filename}...")
        overlapping_gene_names = read_gene_names_from_intervals_file(overlap_interval_filename)

    print("Input files verified.", file=sys.stderr)

    # Create our discovered transcriptome gtf field map:
    gtf_field_dict = get_gtf_field_val_dict(gtf_file)

    # Create our gencode transcriptome gtf field map:
    gencode_gtf_field_dict = get_gtf_field_val_dict(gencode_reference_gtf)

    # Let's read in the equivalence class defs here:
    tx_eq_class_def_map = parse_tx_eq_class_defs_file(tx_eq_class_definitions)
    gene_eq_class_def_map = parse_gene_eq_class_assignments_file(gene_eq_class_assignments)

    # Now read the full eq class associations:
    read_eq_class_map = parse_tx_eq_class_assignments(tx_eq_class_assignments)

    # Create our anndata objects from the given data:
    print("Creating master anndata objects from transcripts counts data...", file=sys.stderr)
    master_adata = create_combined_anndata(
        input_tsv, tx_eq_class_def_map, gene_eq_class_def_map, read_eq_class_map,
        gtf_field_dict, gencode_gtf_field_dict, overlapping_gene_names, overlap_intervals_label
    )

    # Write our data out as pickles:
    print("Pickling data...", file=sys.stderr)
    pickle.dump(master_adata, open(f"{out_prefix}_tx_gene_counts_adata.pickle", "wb"))

    # Write our data as h5ad files:
    print("Writing data to h5ad file...", file=sys.stderr)
    master_adata.write(f"{out_prefix}_tx_gene_counts_adata.h5ad")

    print("Done!", file=sys.stderr)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=f"Creates anndata objects from the given count TSV and a GTF file.",
        epilog="The input TSV should have been created by create_count_matrix_from_annotated_bam.py and therefore"
               "should be of the form:"
               "Gene/TX    CBC    UMI    Count"
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-t', '--tsv',
                               help='TSV file containing gene/transcript counts.',
                               required=True)
    requiredNamed.add_argument('-g', '--gtf',
                               help='GTF file containing gene annotations.  Can be from either Gencode or stringtie.',
                               required=True)
    requiredNamed.add_argument('--tx-eq-class-definitions',
                               help='Transcript equivalence class definitions TSV file.',
                               required=True)
    requiredNamed.add_argument('--tx-eq-class-assignments',
                               help='Transcript equivalence class assignments TSV file.',
                               required=True)
    requiredNamed.add_argument('--gene-eq-class-definitions',
                               help='Gene equivalence class definitions TSV file.',
                               required=True)
    requiredNamed.add_argument('--gene-eq-class-assignments',
                               help='Gene equivalence class assignments TSV file.',
                               required=True)
    requiredNamed.add_argument('-o', '--out-base-name',
                               help='Base name for the output files',
                               required=True)

    parser.add_argument("--overlap-intervals", 
                        help="TSV/Interval list file containing intervals on which to mark transcripts as overlapping.",
                        type=str)
    parser.add_argument("--overlap-interval-label", 
                        help="Label to add to all overlapping intervals.", 
                        type=str, 
                        default="")

    parser.add_argument("--gencode-reference-gtf",
                        help="Gencode GTF file to use to disambiguate the annotations in the given gtf file.",
                        type=str)

    args = parser.parse_args()
    main(args.tsv, args.gtf, args.out_base_name,
         args.tx_eq_class_definitions,
         args.tx_eq_class_assignments,
         args.gene_eq_class_definitions,
         args.gene_eq_class_assignments,
         args.overlap_intervals, args.overlap_interval_label,
         args.gencode_reference_gtf)
