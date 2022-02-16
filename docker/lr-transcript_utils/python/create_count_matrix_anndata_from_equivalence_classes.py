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
                            gtf_field_dict, overlap_intervals=None,
                            overlap_intervals_label="overlaps_intervals_of_interest",
                            gencode_reference_gtf=None,
                            force_recount=False,
                            force_overwrite_gencode_overlaps=False):

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

    cell_barcode_to_tx_to_umi_dict = dict()
    cell_barcode_to_tx_count_dict = dict()

    pickle_file_name = os.path.splitext(os.path.basename(input_tsv))[0] + ".tx_raw_count_matrix.pickle"

    if not force_recount and os.path.exists(pickle_file_name):
        print(f"Loading count map from {pickle_file_name}...", end="\t", file=sys.stderr)
        cell_barcode_to_tx_count_dict = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        # Get our cell tx counts:
        with open(input_tsv, "r") as f, tqdm(desc="Processing Raw Cell Counts", unit=" count") as pbar:
            tsv_file = csv.reader(f, delimiter="\t")
            next(tsv_file)
            for row in tsv_file:
                tx_col = row[0]
                pipe_pos = tx_col.find("|")
                if pipe_pos >= 0:
                    tx_id = tx_col[:tx_col.find("|")]
                else:
                    tx_id = tx_col

                cell_barcode = row[1]
                umi = row[2]

                # Handle the Transcript names:
                if cell_barcode not in cell_barcode_to_tx_to_umi_dict:
                    cell_barcode_to_tx_to_umi_dict[cell_barcode] = dict()
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id] = set()
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id].add(umi)

                    cell_barcode_to_tx_count_dict[cell_barcode] = dict()
                    cell_barcode_to_tx_count_dict[cell_barcode][tx_id] = 1

                elif tx_id not in cell_barcode_to_tx_to_umi_dict[cell_barcode]:
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id] = set()
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id].add(umi)

                    cell_barcode_to_tx_count_dict[cell_barcode][tx_id] = 1

                elif umi not in cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id]:
                    cell_barcode_to_tx_to_umi_dict[cell_barcode][tx_id].add(umi)

                    cell_barcode_to_tx_count_dict[cell_barcode][tx_id] += 1

                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(cell_barcode_to_tx_count_dict, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    ##################################################################################################################
    # Write our cell TX counts to our adata object:

# __   __               _               _   _
# \ \ / /__  _   _     / \   _ __ ___  | | | | ___ _ __ ___
#  \ V / _ \| | | |   / _ \ | '__/ _ \ | |_| |/ _ \ '__/ _ \
#   | | (_) | |_| |  / ___ \| | |  __/ |  _  |  __/ | |  __/
#   |_|\___/ \__,_| /_/   \_\_|  \___| |_| |_|\___|_|  \___|

    # Create unique row / column identifiers into which to aggregate data:
    cell_barcodes = np.array(list(cell_barcode_to_tx_count_dict.keys()))
    tx_ids = np.unique(np.array(list(gtf_field_dict.keys())))

    # Populate the count matrix:
    pickle_file_name = os.path.splitext(os.path.basename(input_tsv))[0] + ".cell_transcript_count_matrix.pickle"
    if not force_recount and os.path.exists(pickle_file_name):
        print(f"Loading count map from {pickle_file_name}...", end="\t", file=sys.stderr)
        count_mat = pickle.load(open(pickle_file_name, "rb"))
        print("Done!", file=sys.stderr)
    else:
        count_mat = scipy.sparse.lil_matrix((len(cell_barcodes), len(tx_ids)), dtype=np.uint32)

        tx_id_index_dict = {name: i for i, name in enumerate(tx_ids)}
        with tqdm(desc=f"Creating cell transcript count matrix", unit=" cell",
                  total=len(cell_barcode_to_tx_count_dict)) as pbar:

            for i, (cb, counts_dict) in enumerate(cell_barcode_to_tx_count_dict.items()):
                # Put the counts for each transcript in the right indices:
                for tx_id, count in counts_dict.items():
                    count_mat[i, tx_id_index_dict[tx_id]] = count
                pbar.update(1)

        print("Pickling data...", file=sys.stderr)
        pickle.dump(count_mat, open(pickle_file_name, "wb"))
        print("Done!", file=sys.stderr)

    ############################################################################################################
    # Now we set up the variables that we're going to apply to each observation:
    # We'll try two different ways - one for Gencode GTFs and one for stringtie GTFs:

    # First we handle our overlapping tests:
    # Mark the transcripts that overlap our intervals:
    if overlap_intervals:
        print(f"Determining overlap intervals for the given interval list... ", end="", file=sys.stderr)
        tx_overlap_flags = [
            interval_overlaps_any_in_interval_list(
                gtf_field_dict[tx_id][CONTIG_FIELD],
                gtf_field_dict[tx_id][START_FIELD],
                gtf_field_dict[tx_id][END_FIELD],
                overlap_intervals
            ) for tx_id in tx_ids
        ]
        print("Done!", file=sys.stderr)

    is_gencode = True
    try:
        # Try Gencode first:
        transcript_names = [gtf_field_dict[tx_id][GENCODE_TX_NAME_FIELD] for tx_id in tx_ids]
        gene_names = [gtf_field_dict[tx_id][GENCODE_GENE_NAME_FIELD] for tx_id in tx_ids]
        gene_ids = [gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]

        de_novo_gene_ids = ["N/A"] * len(transcript_names)
        de_novo_transcript_ids = ["N/A"] * len(transcript_names)
        is_de_novo = [False] * len(transcript_names)
        is_gene_id_ambiguous = [False] * len(transcript_names)

    except KeyError:
        # We must have stringtie data:

        de_novo_transcript_ids = tx_ids
        de_novo_gene_ids = [gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]

        # TX Names for stringtie are the same as TX IDs:
        transcript_names = tx_ids

        # Get our fields for existing transcripts:
        gene_names = [gtf_field_dict[tx_id][STRINGTIE_GENE_NAME_FIELD] if STRINGTIE_GENE_NAME_FIELD in gtf_field_dict[tx_id] else gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]
        gene_ids = [gtf_field_dict[tx_id][STRINGTIE_GENE_ID_FIELD] if STRINGTIE_GENE_ID_FIELD in gtf_field_dict[tx_id] else gtf_field_dict[tx_id][GENE_ID_FIELD] for tx_id in tx_ids]

        # Get our new transcript ids last so we can cue off of them earlier:
        tx_ids = [gtf_field_dict[tx_id][STRINGTIE_TX_ID_FIELD] if STRINGTIE_GENE_NAME_FIELD in gtf_field_dict[tx_id] else tx_id for tx_id in tx_ids]

        # Mark or de novo transcripts:
        is_de_novo = [tx_ids[i] == de_novo_transcript_ids[i] for i in range(len(tx_ids))]

        # Now we have to clean up the gene names.
        # If a stringtie gene overlaps a known gene at any point, we must use the known gene name.
        de_novo_gene_id_to_canonical_gene_name_dict = dict()
        for de_novo_gene_id, canonical_gene_id, canonical_gene_name in zip(de_novo_gene_ids, gene_ids, gene_names):
            if de_novo_gene_id != canonical_gene_id:
                try:
                    de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id].append((canonical_gene_id, canonical_gene_name))
                except KeyError:
                    de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id] = \
                        [(canonical_gene_id, canonical_gene_name)]

        # Make sure we don't have any ambiguous mappings:
        print("Ambiguously / multi-mapped:", file=sys.stderr)
        multimapped = 0
        ambiguous_de_novo_gene_set = set()
        for de_novo_gene_id, gene_info_tuple_list in de_novo_gene_id_to_canonical_gene_name_dict.items():
            if len(gene_info_tuple_list) > 1:
                conflicting = False
                for canonical_gene_id, _ in gene_info_tuple_list:
                    if canonical_gene_id != gene_info_tuple_list[0][0]:
                        conflicting = True
                        break
                if conflicting:
                    multimapped += 1
                    print(f"    {de_novo_gene_id}")
                    gene_info_tuple_count_dict = dict()
                    for canonical_gene_id, canonical_gene_name in gene_info_tuple_list:
                        try:
                            gene_info_tuple_count_dict[canonical_gene_id] += 1
                        except KeyError:
                            gene_info_tuple_count_dict[canonical_gene_id] = 1
                    for canonical_gene_id, canonical_gene_name in set(gene_info_tuple_list):
                        print(f"        {canonical_gene_id} -> {canonical_gene_name} ({gene_info_tuple_count_dict[canonical_gene_id]})")
                    ambiguous_de_novo_gene_set.add(de_novo_gene_id)
        print(f"Num ambiguously / multi-mapped = {multimapped}", file=sys.stderr)

        is_gene_id_ambiguous = [True if dngid in ambiguous_de_novo_gene_set else False for dngid in de_novo_gene_ids]

        # Now let's relabel the genes we know are not ambiguous:
        print("Relabeling unambiguous gene IDs and names...", end="", file=sys.stderr)
        num_relabeled = 0
        for i in range(len(is_de_novo)):
            if is_de_novo[i] and not is_gene_id_ambiguous[i]:
                de_novo_gene_id = de_novo_gene_ids[i]

                fixed_gene_id = None
                fixed_gene_name = None
                try:
                    fixed_gene_id = de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id][0][0]
                    fixed_gene_name = de_novo_gene_id_to_canonical_gene_name_dict[de_novo_gene_id][0][1]
                except KeyError:
                    pass

                if fixed_gene_id and fixed_gene_name:
                    gene_ids[i] = fixed_gene_id
                    gene_names[i] = fixed_gene_name
                    num_relabeled += 1
        print("Done!", file=sys.stderr)
        print(f"Successfully relabeled {num_relabeled} unambiguous gene names / IDs", file=sys.stderr)

        # Let downstream processing know we're not gencode:
        is_gencode = False

    # Create our anndata object now:
    count_adata = anndata.AnnData(count_mat.tocsr())

    # Add our variables:
    col_df = pd.DataFrame()
    col_df["transcript_ids"] = tx_ids
    col_df["gene_ids"] = gene_ids
    col_df["gene_names"] = gene_names
    col_df["transcript_names"] = transcript_names

    col_df["de_novo_gene_ids"] = de_novo_gene_ids
    col_df["de_novo_transcript_ids"] = de_novo_transcript_ids
    col_df["is_de_novo"] = is_de_novo
    col_df["is_gene_id_ambiguous"] = is_gene_id_ambiguous

    # If we're doing interval overlaps, add our label:
    if overlap_intervals:
        print(f"Adding {overlap_intervals_label} column for overlapping intervals...", file=sys.stderr, end="")
        col_df[f"{overlap_intervals_label}"] = tx_overlap_flags
        print("Done!", file=sys.stderr)

    if gencode_reference_gtf:
        if force_overwrite_gencode_overlaps:
            print(f"Adding gencode overlapping gene names (FORCED)...", file=sys.stderr, end="")
            col_df["gencode_overlap_gene_names"] = gene_names
            col_df["gencode_overlap_gene_ids"] = gene_ids
            col_df["is_gencode_gene_overlap_ambiguous"] = [False] * len(gene_names)
            print("Done!", file=sys.stderr)
        else:
            print(f"Adding gencode overlapping gene names...", file=sys.stderr)
            gencode_field_val_dict = get_gtf_field_val_dict(gencode_reference_gtf, entry_type_filter=GENE_ENTRY_STRING)
            print(f"Assigning overlapping genes...", file=sys.stderr)
            raw_overlap_gene_names, raw_overlap_gene_ids, raw_ambiguity_markers = get_approximate_gencode_gene_assignments(gtf_field_dict, gencode_field_val_dict)

            # Reorder by tx_id:
            overlap_gene_names = [None] * len(raw_overlap_gene_names)
            overlap_gene_ids = [None] * len(raw_overlap_gene_ids)
            ambiguity_markers = np.empty(len(raw_ambiguity_markers), dtype=bool)
            for i, tx in enumerate(gtf_field_dict.keys()):
                indx = np.where(de_novo_transcript_ids == tx)[0]
                if len(indx) > 1:
                    raise RuntimeError(f"Error: transcript appears more than once: {tx} ({i}): {indx}")
                # TODO: Remove debugging code:
                # elif len(indx) == 0:
                #     print(f"No TX assignment for {tx} ({i})", file=sys.stderr)
                #     continue

                indx = indx[0]
                print(f"TX Assignment: {tx} ({i}): {indx} - {raw_overlap_gene_names[i]} <{raw_overlap_gene_ids[i]}>", file=sys.stderr)
                overlap_gene_names[indx] = raw_overlap_gene_names[i]
                overlap_gene_ids[indx] = raw_overlap_gene_ids[i]
                ambiguity_markers[indx] = raw_ambiguity_markers[i]

            # TODO: Remove debugging code:
            # for i, v in enumerate(overlap_gene_names):
            #     if v is None:
            #         overlap_gene_names[i] = "None"
            #
            # for i, v in enumerate(overlap_gene_ids):
            #     if v is None:
            #         overlap_gene_ids[i] = "None"

            col_df["gencode_overlap_gene_names"] = overlap_gene_names
            col_df["gencode_overlap_gene_ids"] = overlap_gene_ids
            col_df["is_gencode_gene_overlap_ambiguous"] = ambiguity_markers
            print(f"Num overlap assignments: {len([n for n in overlap_gene_names if (n is not None) and not n.startswith('STRG')])}", file=sys.stderr)
            print(f"Num ambiguous assignments: {len(np.where(ambiguity_markers == True)[0])}", file=sys.stderr)
            print("Done!", file=sys.stderr)

    # Assign the data to our anndata object:
    count_adata.var = col_df

    if is_gencode:
        count_adata.var_names = transcript_names
    else:
        count_adata.var_names = tx_ids

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


def parse_eq_class_defs_file(file_name):
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


def parse_tx_eq_class_assignments(file_name):
    read_eq_class_map = dict()
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            read_name, tx_eq_class, gene_assignment = line.strip().split("\t")
            read_eq_class_map[read_name] = tuple(tx_eq_class, gene_assignment)

    return read_eq_class_map


def main(input_tsv, gtf_file, out_prefix,
         tx_eq_class_definitions,
         tx_eq_class_assignments,
         gene_eq_class_definitions,
         gene_eq_class_assignments,
         overlap_interval_filename=None, overlap_intervals_label=None,
         gencode_reference_gtf=None, force_overwrite_gencode_overlaps=False):

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

    overlap_intervals = None
    if overlap_interval_filename:
        print("Verifying contents of {overlap_interval_filename}...")
        overlap_intervals = read_intervals_from_tsv(overlap_interval_filename)

    print("Input files verified.", file=sys.stderr)

    # Create our gtf field map:
    gtf_field_dict = get_gtf_field_val_dict(gtf_file)

    # Let's read in the equivalence class defs here:
    tx_eq_class_def_map = parse_eq_class_defs_file(tx_eq_class_definitions)
    gene_eq_class_def_map = parse_eq_class_defs_file(gene_eq_class_definitions)

    # Now read the full eq class associations:
    read_eq_class_map = parse_tx_eq_class_assignments(tx_eq_class_assignments)

    # Create our anndata objects from the given data:
    print("Creating master anndata objects from transcripts counts data...", file=sys.stderr)
    master_adata = create_combined_anndata(
        input_tsv, tx_eq_class_def_map, gene_eq_class_def_map, read_eq_class_map,
        gtf_field_dict, overlap_intervals, overlap_intervals_label,
        gencode_reference_gtf, force_overwrite_gencode_overlaps=force_overwrite_gencode_overlaps
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

    parser.add_argument("--force-overwrite-gencode-overlaps",
                        help="Forces the values in `gene_names` and `gene_ids` to overwrite `gencode_overlap_gene_names`"
                             " and `gencode_overlap_gene_ids` respectively.  This is only valid if given with "
                             "--gencode-reference-gtf.",
                        action='store_true')

    args = parser.parse_args()
    main(args.tsv, args.gtf, args.out_base_name,
         args.tx_eq_class_definitions,
         args.tx_eq_class_assignments,
         args.gene_eq_class_definitions,
         args.gene_eq_class_assignments,
         args.overlap_intervals, args.overlap_interval_label,
         args.gencode_reference_gtf, args.force_overwrite_gencode_overlaps)
