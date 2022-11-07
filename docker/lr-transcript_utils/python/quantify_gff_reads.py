#!/usr/bin/env scripts

import time
import sys
import pysam
import argparse

import numpy as np
import networkx as nx

from collections import defaultdict
from tqdm import tqdm

##########################################################################################
##########################################################################################
##########################################################################################

CC_OFFSET = 1

cc_eq  = 0   + CC_OFFSET
cc_c   = 1   + CC_OFFSET
cc_k   = 2   + CC_OFFSET
cc_m   = 3   + CC_OFFSET
cc_n   = 4   + CC_OFFSET
cc_j   = 5   + CC_OFFSET
cc_e   = 6   + CC_OFFSET
cc_o   = 7   + CC_OFFSET
cc_s   = 8   + CC_OFFSET
cc_x   = 9   + CC_OFFSET
cc_i   = 10  + CC_OFFSET
cc_y   = 11  + CC_OFFSET
cc_p   = 12  + CC_OFFSET
cc_r   = 13  + CC_OFFSET
cc_u   = 14  + CC_OFFSET
cc_dot = 15  + CC_OFFSET

class_code_lookup = {
    "=" : cc_eq,
    "c" : cc_c,
    "k" : cc_k,
    "m" : cc_m,
    "n" : cc_n,
    "j" : cc_j,
    "e" : cc_e,
    "o" : cc_o,
    "s" : cc_s,
    "x" : cc_x,
    "i" : cc_i,
    "y" : cc_y,
    "p" : cc_p,
    "r" : cc_r,
    "u" : cc_u,
    "." : cc_dot,
}

class_code_num_lookup = {
   cc_eq  :  "=",
   cc_c   :  "c",
   cc_k   :  "k",
   cc_m   :  "m",
   cc_n   :  "n",
   cc_j   :  "j",
   cc_e   :  "e",
   cc_o   :  "o",
   cc_s   :  "s",
   cc_x   :  "x",
   cc_i   :  "i",
   cc_y   :  "y",
   cc_p   :  "p",
   cc_r   :  "r",
   cc_u   :  "u",
   cc_dot :  ".",
}

gencode_node_filter = lambda n: n.startswith('ENST')
stringtie2_node_filter = lambda n: n.startswith('STRG')
mas_seq_node_filter = lambda n: n.startswith('m6402')

gencode_gene_name_filter = lambda n: not n.startswith("STRG") and not n.startswith("m6402")

##########################################################################################
##########################################################################################
##########################################################################################


def get_field_size(values):
    max_field = -np.inf
    for f in values:
        if not isinstance(f, list) and not isinstance(f, (list, np.ndarray)):
            if f > max_field:
                max_field = f
    nd = int(np.ceil(np.log10(max_field)))
    return nd


def print_stats(a, label=None):
    prec = ".4f"
    int_types = {
        int,
        np.int,
        np.intc,
        np.intp,
        np.int0,
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        np.integer,
        float,
        np.uint,
        np.uintc,
        np.uintp,
        np.uint0,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.unsignedinteger
    }

    quantiles_to_calc = [0.1, 0.25, 0.75, 0.9]
    quantiles = np.quantile(a, quantiles_to_calc)

    # Construct our data dictionary in a way that will display our stats in
    # a sensical manner:
    data = {
        "#": len(a),
        "Min": np.min(a),
        "Max": np.max(a),
    }
    for i in range(int(len(quantiles_to_calc) / 2)):
        q = quantiles_to_calc[i]
        v = quantiles[i]
        data[f"Q{int(q * 100)}"] = v

    data["Mean"] = np.mean(a)
    data["Median"] = np.median(a)

    for i in range(int(len(quantiles_to_calc) / 2), len(quantiles_to_calc)):
        q = quantiles_to_calc[i]
        v = quantiles[i]
        data[f"Q{int(q * 100)}"] = v

    data["Stdev"] = np.std(a)

    nd = get_field_size(data.values())

    field_width = max([len(k) for k in data.keys()])

    if label:
        print(f"{label}:")
        print("-" * (len(label) + 1))

    for k, v in data.items():
        if type(v) in int_types:
            print(f"{k:{field_width}}: {v:{nd}}")
        else:
            if v == 0:
                nfd = 1
            else:
                nfd = int(np.ceil(np.log10(v)))
            spacer = " " * (nd - nfd)
            print(f"{k:{field_width}}: {spacer}{v:{prec}}")


def _line_count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)


def count_lines_in_file(file_name):
    # Get quick count of lines in the file:
    with open(file_name, 'rb') as fp:
        c_generator = _line_count_generator(fp.raw.read)
        num_lines = sum(buffer.count(b'\n') for buffer in c_generator) - 1
    return num_lines


def has_edge(nx_graph, start_node, end_node, class_code):
    edge_exists = False
    for src, dest, key, dat in nx_graph.out_edges([start_node], data=True, keys=True):
        if src == start_node and dest == end_node and dat["class_code"] == class_code:
            edge_exists = True
            break

    return edge_exists


def add_refmap_to_graph(refmap_file, nx_graph=None, gene_name_gene_id_map=None, do_eq_filter=False):
    if not nx_graph:
        nx_graph = nx.MultiDiGraph()

    num_lines = count_lines_in_file(refmap_file)

    with open(refmap_file, 'r') as f:
        with tqdm(desc=f"Parsing refmap for graph links...", unit=" line", total=num_lines) as pbar:
            for i, line in enumerate(f):
                if i == 0:
                    continue

                ref_gene, ref_tx, cc, query_txs = line.strip().split('\t')

                ref_gene = ref_gene.strip()
                ref_tx = ref_tx.strip()
                cc = cc.strip()

                if do_eq_filter and cc != "=":
                    continue

                if gene_name_gene_id_map and ref_gene in gene_name_gene_id_map:
                    ref_gene = gene_name_gene_id_map[ref_gene]

                if ref_tx not in nx_graph:
                    nx_graph.add_node(ref_tx, gene_name_set={ref_gene})

                nx_graph.nodes[ref_tx]["gene_name_set"].add(ref_gene)

                for q in query_txs.split(","):
                    qg, qtx = q.split("|")
                    qg = qg.strip()
                    qtx = qtx.strip()

                    if gene_name_gene_id_map and qg in gene_name_gene_id_map:
                        qg = gene_name_gene_id_map[qg]

                    if qtx not in nx_graph:
                        nx_graph.add_node(qtx, gene_name_set={qg})

                    nx_graph.nodes[qtx]["gene_name_set"].add(qg)

                    # Add an edge from qtx to ref_tx if it doesn't already exist:
                    if not nx_graph.has_edge(qtx, ref_tx, key=cc):
                        nx_graph.add_edge(qtx, ref_tx, key=cc, class_code=cc)

                pbar.update(1)
    return nx_graph


def add_tmap_to_graph(tmap_file, nx_graph=None, gene_name_gene_id_map=None, add_link_to_major_isoform_id=False):
    if not nx_graph:
        nx_graph = nx.MultiDiGraph()

    num_lines = count_lines_in_file(tmap_file)

    major_isoform_map = dict()

    with open(tmap_file, 'r') as f:
        with tqdm(desc=f"Parsing tmap for graph links...", unit=" line", total=num_lines) as pbar:
            for i, line in enumerate(f):
                if i == 0:
                    continue

                ref_gene, ref_tx, cc, query_gene, query_tx, num_exons, fpkm, tpm, cov, length, major_isoform_id, ref_match_len = line.strip().split(
                    '\t')

                ref_gene = ref_gene.strip()
                ref_tx = ref_tx.strip()
                cc = cc.strip()
                query_gene = query_gene.strip()
                query_tx = query_tx.strip()
                num_exons = num_exons.strip()
                fpkm = fpkm.strip()
                tpm = tpm.strip()
                cov = cov.strip()
                length = length.strip()
                major_isoform_id = major_isoform_id.strip()
                ref_match_len = ref_match_len.strip()

                #############################################
                # Handle the reference gene / transcript:
                if gene_name_gene_id_map and ref_gene in gene_name_gene_id_map:
                    ref_gene = gene_name_gene_id_map[ref_gene]

                if ref_tx not in nx_graph:
                    nx_graph.add_node(ref_tx, gene_name_set={ref_gene})

                nx_graph.nodes[ref_tx]["gene_name_set"].add(ref_gene)

                #############################################
                # Handle the query gene / transcript:
                if gene_name_gene_id_map and query_gene in gene_name_gene_id_map:
                    query_gene = gene_name_gene_id_map[query_gene]

                if query_tx not in nx_graph:
                    if not mas_seq_node_filter(query_gene):
                        nx_graph.add_node(query_tx, gene_name_set={query_gene})
                    else:
                        nx_graph.add_node(query_tx, gene_name_set=set())

                if not mas_seq_node_filter(query_gene):
                    nx_graph.nodes[query_tx]["gene_name_set"].add(query_gene)

                # Add an edge from query_tx to ref_tx:
                if not nx_graph.has_edge(query_tx, ref_tx, key=cc):
                    nx_graph.add_edge(query_tx, ref_tx, key=cc, class_code=cc)

                #############################################
                # track our major isoform ids if we need to:
                if add_link_to_major_isoform_id:
                    major_isoform_map[query_tx] = major_isoform_id

                pbar.update(1)

    if add_link_to_major_isoform_id:
        # OK, we need to add in our links now that we have all the IDs:
        with tqdm(desc=f"Adding links to major isoform IDs...", unit=" query txs",
                  total=len(major_isoform_map)) as pbar:
            for query_tx, major_isoform_id in major_isoform_map.items():
                if not nx_graph.has_edge(query_tx, major_isoform_id, key="="):
                    nx_graph.add_edge(query_tx, major_isoform_id, key="=")
                pbar.update(1)

    return nx_graph


def get_gencode_gene_id_gene_name_map(gencode_file):
    num_lines = count_lines_in_file(gencode_file)

    gene_name_gene_id_map = dict()
    gene_id_gene_name_map = dict()

    with open(gencode_file, 'r') as f:
        with tqdm(desc=f"Getting gencode gene ID->gene name.", unit=" line", total=num_lines) as pbar:
            for i, line in enumerate(f):
                if line.startswith("#"):
                    pbar.update(1)
                    continue

                contig, source, feature_type, start, end, _, strand, _, extra_fields = line.strip().split('\t')

                if feature_type != "gene":
                    pbar.update(1)
                    continue

                # We know we're in a gene.  Now we get the gene id and gene name:
                gene_name = None
                gene_id = None
                for field in extra_fields.split(";"):
                    field_name, field_value = field.strip().split(" ")
                    if field_name == "gene_id":
                        gene_id = field_value.replace('"', "")
                    elif field_name == "gene_name":
                        gene_name = field_value.replace('"', "")

                    if gene_name and gene_id:
                        break

                gene_name_gene_id_map[gene_name] = gene_id
                gene_id_gene_name_map[gene_id] = gene_name

                pbar.update(1)

    return gene_name_gene_id_map, gene_id_gene_name_map


def get_mas_seq_read_tx_assignments(graph):
    # Getting Gencode TX Assignments:
    # 1) Find all single "=" assignments (i.e. a string of '=' to an 'ENST' resolution)
    # 2) Everything else:
    #      Collect all genes, transcripts, and the class codes and make a list of them.

    num_mas_seq_nodes = nx.subgraph_view(graph, filter_node=mas_seq_node_filter).number_of_nodes()

    # 1) Get all definitive assignments:
    mas_seq_tx_equivalence_classes = dict()
    with tqdm(desc=f"Extracting exact mas-seq -> gencode tx", unit=" node", total=num_mas_seq_nodes) as pbar:
        for mas_seq_node in nx.subgraph_view(graph, filter_node=mas_seq_node_filter).nodes:

            # Get all '=' nodes connected to this one:
            node_traversal_list = [(mas_seq_node, [])]
            nodes_seen = dict()
            while len(node_traversal_list) > 0:
                n, path = node_traversal_list.pop()

                if n in nodes_seen:
                    continue

                # TODO: We don't need to prefer `=` edges anymore!  The EQ Class is a subgraph, so to be consistent we should not do this.
                for src, dest, d in graph.out_edges([n], data=True):
                    cc = d["class_code"]
                    if cc == "=":
                        new_path = [p for p in path]
                        new_path.append(cc)
                        node_traversal_list.append((dest, new_path))

                nodes_seen[n] = path

            # Now we have a list of nodes directly connected to our mas-seq node.
            # Pull out all tx definitions from them into a set.

            tx_assignments = set()
            for n, path in nodes_seen.items():
                # TODO: Shouldn't this be gencode OR stringtie nodes?
                if gencode_node_filter(n):
                    tx_assignments.add((n, '='))

            # Only make an assignment here if we have one:
            if len(tx_assignments) == 1:
                mas_seq_tx_equivalence_classes[mas_seq_node] = tx_assignments
            pbar.update(1)

    # 2) Get all other tx info:
    num_inexact = num_mas_seq_nodes - len(mas_seq_tx_equivalence_classes)
    with tqdm(desc=f"Extracting in-exact mas-seq -> gencode tx", unit=" node", total=num_inexact) as pbar:
        for mas_seq_node in nx.subgraph_view(graph,
                                             filter_node=lambda n: mas_seq_node_filter(n) and
                                                                   n not in mas_seq_tx_equivalence_classes).nodes:

            # Get all nodes connected to this one on out edges
            # for each of them, record the edge class code and destination:
            tx_assignments = set()
            for src, dest, d in graph.out_edges([mas_seq_node], data=True):
                cc = d["class_code"]
                tx_assignments.add((dest, cc))

            # Check here for null (-) assignments.
            # If we have exactly 1 null assignment and 1 of any other kind of assignment
            # we ignore the null assignment.
            #
            # This resolves the case where a read maps to the genome but not to stringtie and visa-versa.
            # For example:
            #
            #     Gencode / Read Tmap:
            #          ENSG00000000460.17  ENST00000498289.5 i m64020e_210505_070001/453446566/ccs m64020e_210505_070001/453446566/ccs 1 0.000000  0.000000  0.000000  777 m64020e_210505_070001/453446566/ccs 3849
            #     ST2 / Read Tmap:
            #          - - u m64020e_210505_070001/453446566/ccs m64020e_210505_070001/453446566/ccs 1 0.000000  0.000000  0.000000  777 m64020e_210505_070001/453446566/ccs -
            #
            # TODO: This should actually drop the null assignments in the case where there are ANY other assignments.
            num_null = 0
            for dest, cc in tx_assignments:
                if dest == "-":
                    num_null += 1
            if num_null == 1 and len(tx_assignments) == 2:
                new_tx_assignments = set()
                for dest, cc in tx_assignments:
                    if dest != "-":
                        new_tx_assignments.add((dest, cc))
                tx_assignments = new_tx_assignments

            mas_seq_tx_equivalence_classes[mas_seq_node] = tx_assignments
            pbar.update(1)

    return mas_seq_tx_equivalence_classes


def get_mas_seq_read_gene_assignments(graph):
    num_mas_seq_nodes = nx.subgraph_view(graph, filter_node=mas_seq_node_filter).number_of_nodes()

    mas_seq_to_gencode_gene = dict()
    with tqdm(desc=f"Extracting mas-seq gene map", unit=" node", total=num_mas_seq_nodes) as pbar:
        for mas_seq_node in nx.subgraph_view(graph, filter_node=mas_seq_node_filter).nodes:

            # Get all nodes connected to this one:
            node_traversal_list = [(mas_seq_node, [])]
            nodes_seen = dict()
            while len(node_traversal_list) > 0:
                n, path = node_traversal_list.pop()

                if n in nodes_seen:
                    continue

                # Only continue on if we aren't at a gencode node:
                if not gencode_node_filter(n):
                    for src, dest, d in graph.out_edges([n], data=True):
                        if dest not in nodes_seen:
                            new_path = [p for p in path]
                            new_path.append(d["class_code"])
                            node_traversal_list.append((dest, new_path))

                nodes_seen[n] = path

            # Now we have a list of nodes directly connected to our mas-seq node.
            # Pull out all 'ENSG' gene definitions from them into a set.

            paths_have_eq = False
            connected_gencode_genes = defaultdict(list)
            for n, path in nodes_seen.items():
                for gene_name in graph.nodes[n]["gene_name_set"]:
                    if gencode_gene_name_filter(gene_name):
                        if "=" in path:
                            paths_have_eq = True
                        connected_gencode_genes[gene_name].append(path)

            # Now we can refine our entries by looking at the paths.
            # 1) If there are no gene names, we're done.

            # 2) If there is a single gene name, we're done.
            # 3) If there are more than 1 gene name:
            #    1) If there are paths with `=` in them, we remove all other paths.
            #    2) Keep all remaining gene names.

            gene_assignments = set()
            for gencode_gene_name, path_list in connected_gencode_genes.items():
                if paths_have_eq:
                    is_eq_path = False
                    for path in path_list:
                        if "=" in path:
                            is_eq_path = True
                            break
                    if not is_eq_path:
                        continue
                gene_assignments.add(gencode_gene_name)

            # Check here for null (-) assignments.
            # If we have exactly 1 null assignment and 1 of any other kind of assignment
            # we ignore the null assignment.
            #
            # This resolves the case where a read maps to the genome but not to stringtie and visa-versa.
            # For example:
            #
            #     Gencode / Read Tmap:
            #          ENSG00000000460.17  ENST00000498289.5 i m64020e_210505_070001/453446566/ccs m64020e_210505_070001/453446566/ccs 1 0.000000  0.000000  0.000000  777 m64020e_210505_070001/453446566/ccs 3849
            #     ST2 / Read Tmap:
            #          - - u m64020e_210505_070001/453446566/ccs m64020e_210505_070001/453446566/ccs 1 0.000000  0.000000  0.000000  777 m64020e_210505_070001/453446566/ccs -
            #
            # TODO: This should actually drop the null assignments in the case where there are ANY other assignments.
            num_null = 0
            for g in gene_assignments:
                if g == "-":
                    num_null += 1
            if num_null == 1 and len(gene_assignments) == 2:
                gene_assignments.remove("-")

            mas_seq_to_gencode_gene[mas_seq_node] = gene_assignments
            pbar.update(1)

    return mas_seq_to_gencode_gene


def update_stringtie2_gene_names(graph):
    # Now resolve stringtie 2 gene names to gencode gene names if they aren't all labeled properly.
    # 1) Get list of stringtie2 genes and their mapping to gencode genes.
    # 2) Remove multi-mapped / ambiguous genes from the list.
    # 3) For each st2 tx:
    #    1) get the st2 gene name
    #    2) Get the gencode mapping for that st2 gene name
    #    3) Get all out edges from the st2 tx node
    #       1) if this st2 node does not have out edges that connect to gencode nodes
    #          1)  Add the gencode name to the st2 node gene name list

    num_st2_nodes = nx.subgraph_view(graph, filter_node=stringtie2_node_filter).number_of_nodes()

    st2_gene_to_gencode_genes = dict()
    with tqdm(desc=f"Extracting raw st2 -> gencode gene", unit=" node", total=num_st2_nodes) as pbar:
        for st2_node in nx.subgraph_view(graph, filter_node=stringtie2_node_filter).nodes:

            gencode_conn_list = []
            has_eq = False
            # Each ST2 node should have out edges directly to a gencode gene:
            for src, dest, d in graph.out_edges([st2_node], data=True):
                if gencode_node_filter(dest):
                    gencode_gene_names = [name for name in graph.nodes[dest]["gene_name_set"] if
                                          gencode_gene_name_filter(name)]
                    for ggn in gencode_gene_names:
                        gencode_conn_list.append((ggn, src, d['class_code']))
                    if d['class_code'] == "=":
                        has_eq = True

            if has_eq:
                gencode_conn_list = [(ggn, src, cc) for ggn, src, cc in gencode_conn_list if cc == "="]

            # Really, this should only ever have 1 name in it:
            for st2_gene_name in graph.nodes[st2_node]["gene_name_set"]:
                if st2_gene_name in st2_gene_to_gencode_genes:
                    st2_gene_to_gencode_genes[st2_gene_name].extend(gencode_conn_list)
                else:
                    st2_gene_to_gencode_genes[st2_gene_name] = gencode_conn_list

            pbar.update(1)

    # Now we have to filter the assignments by the relationships we've accumulated:
    final_st2_gene_to_gencode_genes = dict()
    with tqdm(desc=f"Refining st2 -> gencode gene map", unit=" node", total=len(st2_gene_to_gencode_genes)) as pbar:
        for st2_gene_name, gencode_conn_list in st2_gene_to_gencode_genes.items():

            #TODO: This is probably unnecessary since we're filtering out everything else in the case of an `=` relationship.
            has_eq = any([cc == '=' for ggn, src, cc in gencode_conn_list])

            unique_gencode_genes = set()
            new_conn_list = []
            for ggn, src, cc in gencode_conn_list:
                # If we have an equal relationship here and the current
                # link is not an equals, we ignore it:
                if has_eq and cc != '=':
                    continue

                # If we have already seen this gene name, we can ignore it:
                if ggn in unique_gencode_genes:
                    continue

                # OK, we have something new, we should track it:
                new_conn_list.append((ggn, src, cc))
                unique_gencode_genes.add(ggn)

            final_st2_gene_to_gencode_genes[st2_gene_name] = new_conn_list

            pbar.update(1)

    # Now we add in the new gene name assignments to the ST2 nodes if they are unambiguous:
    with tqdm(desc=f"Updating graph st2 -> gencode gene map", unit=" node", total=num_st2_nodes) as pbar:
        for st2_node in nx.subgraph_view(graph, filter_node=stringtie2_node_filter).nodes:

            # There should only ever be one stringtie 2 gene:
            st2_gene_name = [gn for gn in graph.nodes[st2_node]["gene_name_set"] if gn.startswith("STRG")][0]

            # We know we're working with unique gencode genes:
            if len(final_st2_gene_to_gencode_genes[st2_gene_name]) > 0 and len(
                    final_st2_gene_to_gencode_genes[st2_gene_name][0]) == 1:
                gencode_gene_assignment = final_st2_gene_to_gencode_genes[st2_gene_name][0][0]
                cc = final_st2_gene_to_gencode_genes[st2_gene_name][0][2]

                # Update the node -> gene map:
                s = graph.nodes[st2_node]["gene_name_set"]
                s.add((gencode_gene_assignment, cc))
                graph.nodes[st2_node]["gene_name_set"] = s

            pbar.update(1)

    return graph


def write_read_to_gene_assignment_files(out_base_name, gene_eq_classes, mas_seq_to_gencode_gene):

    # Write out EQ Class label file:
    with open(f"{out_base_name}.gene_equivalence_class_lookup.tsv", 'w') as f:
        with tqdm(desc=f"Writing Gene EQ class label file", unit=" eq class", total=len(gene_eq_classes)) as pbar:
            f.write("#EQ_Class\tGene_Assignments\n")
            for eq_class, indx in gene_eq_classes.items():
                f.write(f"MASG{indx:011d}\t{','.join(eq_class)}\n")
                pbar.update(1)

    # Write out gene assignments:
    with open(f"{out_base_name}.gene_name_assignments.tsv", 'w') as f:
        with tqdm(desc=f"Writing mas-seq -> gene map", unit=" mas read", total=len(mas_seq_to_gencode_gene)) as pbar:
            for k, v in mas_seq_to_gencode_gene.items():
                # Don't write out gene assignments for reads that don't have them:
                if len(v) == 1:
                    if list(v)[0] == "-":
                        # We take the EQ class for this read because it is a solitary read:
                        eq_class_label = gene_eq_classes[tuple([k])]
                        f.write(f"{k}\tMASG{eq_class_label:011d}\n")
                    else:
                        # No eq class necessary.
                        f.write(f"{k}\t{','.join(sorted(v))}\n")
                elif len(v) > 1:
                    # Need to relabel with eq class:
                    eq_class = tuple(sorted(list(v)))
                    eq_class_label = gene_eq_classes[eq_class]
                    f.write(f"{k}\tMASG{eq_class_label:011d}\n")
                pbar.update(1)


def write_read_to_equivalence_class_files(out_base_file_name,
                                          tx_eq_class_defs,
                                          mas_seq_tx_equivalence_classes,
                                          gene_eq_classes,
                                          mas_seq_to_gencode_gene):
    # Write out EQ Class label file:
    with open(f"{out_base_file_name}.equivalence_class_lookup.tsv", 'w') as f:
        with tqdm(desc=f"Writing EQ class label file", unit=" eq class", total=len(tx_eq_class_defs)) as pbar:
            f.write("#EQ_Class\tTranscript_Assignments\n")
            for eq_class, indx in tx_eq_class_defs.items():
                f.write(f"MAST{indx:011d}\t{','.join([str(tx) + ';' + str(cc) for tx, cc in eq_class])}\n")
                pbar.update(1)

    with open(f"{out_base_file_name}.equivalence_classes.tsv", 'w') as f:
        with tqdm(desc=f"Writing mas-seq -> eq class", unit=" mas read",
                  total=len(mas_seq_tx_equivalence_classes)) as pbar:
            f.write("#Read_Name\tEQ_Class\tAssociated_Gene(s)\n")
            for mas_read, tx_set in mas_seq_tx_equivalence_classes.items():
                gene_assignments = tuple(sorted(list(mas_seq_to_gencode_gene[mas_read])))

                if len(gene_assignments) == 0:
                    gene_eq_class_label = ""
                elif len(gene_assignments) == 1:
                    gene_eq_class_label = gene_assignments[0]
                else:
                    gene_eq_class_label = f"MASG{gene_eq_classes[gene_assignments]:011d}"

                tx_eq_class = tuple(sorted(list(tx_set)))
                tx_eq_class_label = tx_eq_class_defs[tx_eq_class]

                f.write(f"{mas_read}\tMAST{tx_eq_class_label:011d}\t{gene_eq_class_label}\n")
                pbar.update(1)


def print_node_info(graph, node_name):
    print(node_name)
    print("Data:")
    for k,v in graph.nodes[node_name].items():
        print(f"\t{k}: {v}")
    print("Out edges:")
    for src, dest, dat in graph.out_edges([node_name], data=True):
        print(f"\t{src}\t-[{dat['class_code']}]-> {dest}")
    print("In edges:")
    for src, dest, dat in graph.in_edges([node_name], data=True):
        print(f"\t{src}\t-[{dat['class_code']}]-> {dest}")


def print_graph_stats(graph):
    node_print_width = int(np.log10(graph.number_of_nodes()) + 1)
    edge_print_width = int(np.log10(graph.number_of_edges()) + 1)
    print_width = max([node_print_width, edge_print_width])

    print(f"Total Transcripts:\t{graph.number_of_nodes():{print_width}d}")
    print(f"Total Graph Edges:\t{graph.number_of_edges():{print_width}d}")
    print()

    print(
        f"Num gencode TX entries:\t{nx.subgraph_view(graph, filter_node=gencode_node_filter).number_of_nodes():{print_width}d}")
    print(
        f"Num st2 TX entries:\t{nx.subgraph_view(graph, filter_node=stringtie2_node_filter).number_of_nodes():{print_width}d}")
    print(
        f"Num MAS-seq TX entries:\t{nx.subgraph_view(graph, filter_node=mas_seq_node_filter).number_of_nodes():{print_width}d}")

    gencode_gene_to_tx_dict = defaultdict(list)
    for gencode_node in nx.subgraph_view(graph, filter_node=gencode_node_filter).nodes:
        for gn in graph.nodes[gencode_node]["gene_name_set"]:
            gencode_gene_to_tx_dict[gn].append(gencode_node)

    st2_gene_to_tx_dict = defaultdict(list)
    for st2_node in nx.subgraph_view(graph, filter_node=stringtie2_node_filter).nodes:
        for gn in graph.nodes[st2_node]["gene_name_set"]:
            st2_gene_to_tx_dict[gn].append(st2_node)

    print()
    print(f"Number of Gencode Genes: {len(gencode_gene_to_tx_dict)}")
    print(f"Number of ST2 Genes: {len(st2_gene_to_tx_dict)}")

    print()
    print("Edge type breakdown:")
    edge_type_dict = defaultdict(int)
    for source, dest, attrs in graph.edges(data=True):
        edge_type_dict[attrs["class_code"]] += 1

    for k, v in edge_type_dict.items():
        print(f"\t{k} ({class_code_lookup[k]})\t{v}")


def main(gencode_gtf,
         st2_gencode_refmap, st2_mas_seq_refmap, gencode_st2_refmap, gencode_mas_seq_refmap,
         st2_gencode_tmap, st2_mas_seq_tmap, gencode_st2_tmap, gencode_mas_seq_tmap,
         out_base_name):

    # Get the mapping from gencode gene name to gene ID so we can clean up the gene names in teh graph:
    gene_name_gene_id_map, gene_id_gene_name_map = get_gencode_gene_id_gene_name_map(gencode_gtf)

    # Create our graph:
    graph = nx.MultiDiGraph()

    # Add refmap files:
    graph = add_refmap_to_graph(st2_gencode_refmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)
    graph = add_refmap_to_graph(st2_mas_seq_refmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)
    graph = add_refmap_to_graph(gencode_st2_refmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)
    graph = add_refmap_to_graph(gencode_mas_seq_refmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)

    # Add tmap files:
    graph = add_tmap_to_graph(st2_gencode_tmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)
    graph = add_tmap_to_graph(st2_mas_seq_tmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)
    graph = add_tmap_to_graph(gencode_st2_tmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)
    graph = add_tmap_to_graph(gencode_mas_seq_tmap, graph, gene_name_gene_id_map=gene_name_gene_id_map)

    # Update StringTie2 gene names with resolved gencode names:
    graph = update_stringtie2_gene_names(graph)

    # Get our gene name assignments:
    mas_seq_to_gencode_gene = get_mas_seq_read_gene_assignments(graph)

    # Get our transcript / equivalence class assignments:
    mas_seq_tx_equivalence_classes = get_mas_seq_read_tx_assignments(graph)

    # Print some stats about our gene and transcript assignments:
    num_mas_seq_nodes = nx.subgraph_view(graph, filter_node=mas_seq_node_filter).number_of_nodes()

    print()
    print(f"Num MAS-seq nodes: {num_mas_seq_nodes}")
    print()

    if num_mas_seq_nodes == 0:
        print("No mas-seq nodes in this dataset.  This may be fine if we're on an alt contig.")
    else:
        print("Gene Stats:")
        mas_gencode_gene_assignment_counts = [len(s) for s in mas_seq_to_gencode_gene.values()]
        print(
            f"Num uniquely assigned MAS-seq reads: {sum([c == 1 for c in mas_gencode_gene_assignment_counts])}/{num_mas_seq_nodes} ({100 * sum([c == 1 for c in mas_gencode_gene_assignment_counts]) / num_mas_seq_nodes:2.4f}%)")
        print(
            f"Num ambiguous/unassigned MAS-seq reads: {sum([c == 0 for c in mas_gencode_gene_assignment_counts])}/{num_mas_seq_nodes} ({100 * sum([c == 0 for c in mas_gencode_gene_assignment_counts]) / num_mas_seq_nodes:2.4f}%)")
        print(
            f"Num multi-assigned MAS-seq reads: {sum([c > 1 for c in mas_gencode_gene_assignment_counts])}/{num_mas_seq_nodes} ({100 * sum([c > 1 for c in mas_gencode_gene_assignment_counts]) / num_mas_seq_nodes:2.4f}%)")
        print()
        print_stats(mas_gencode_gene_assignment_counts, "MAS-seq gene assignment stats")
        print()

        print()

        print("Transcript Stats:")
        mas_gencode_tx_assignment_counts = [len(s) for s in mas_seq_tx_equivalence_classes.values()]
        print(
            f"Num uniquely assigned MAS-seq reads: {sum([c == 1 for c in mas_gencode_tx_assignment_counts])}/{num_mas_seq_nodes} ({100 * sum([c == 1 for c in mas_gencode_tx_assignment_counts]) / num_mas_seq_nodes:2.4f}%)")
        print(
            f"Num ambiguous/unassigned MAS-seq reads: {sum([c == 0 for c in mas_gencode_tx_assignment_counts])}/{num_mas_seq_nodes} ({100 * sum([c == 0 for c in mas_gencode_tx_assignment_counts]) / num_mas_seq_nodes:2.4f}%)")
        print(
            f"Num multi-assigned MAS-seq reads: {sum([c > 1 for c in mas_gencode_tx_assignment_counts])}/{num_mas_seq_nodes} ({100 * sum([c > 1 for c in mas_gencode_tx_assignment_counts]) / num_mas_seq_nodes:2.4f}%)")
        print()
        print_stats(mas_gencode_tx_assignment_counts, "MAS-seq transcript assignment stats")
        print()

    # Calculate unique gene equivalence classes:
    gene_eq_classes = set()
    for mas_read, gene_set in mas_seq_to_gencode_gene.items():
        # For solitary null ('-') gene assignments, we create a separate eq class:
        if len(gene_set) == 1 and list(gene_set)[0] == "-":
            gene_eq_classes.add(tuple([mas_read]))
        elif len(gene_set) > 1:
            gene_eq_classes.add(tuple(sorted(list(gene_set))))
    gene_eq_classes = {eq_class: indx for indx, eq_class in enumerate(tuple(sorted(list(gene_eq_classes))))}

    # Calculate unique transcript equivalence classes:
    with tqdm(desc=f"Creating EQ class labels", unit=" eq class", total=len(mas_seq_tx_equivalence_classes)) as pbar:
        tx_eq_classes = set()
        for eq_class in mas_seq_tx_equivalence_classes.values():
            tx_eq_classes.add(tuple(sorted(list(eq_class))))
            pbar.update(1)
        tx_eq_classes = {eq_class: indx for indx, eq_class in enumerate(tuple(sorted(list(tx_eq_classes))))}

    print(f"Num unique equivalence classes: {len(tx_eq_classes)}")

    # Write our gene assignments:
    write_read_to_gene_assignment_files(out_base_name,
                                        gene_eq_classes,
                                        mas_seq_to_gencode_gene)

    # Write our equivalence classes:
    write_read_to_equivalence_class_files(out_base_name,
                                          tx_eq_classes,
                                          mas_seq_tx_equivalence_classes,
                                          gene_eq_classes,
                                          mas_seq_to_gencode_gene)

    # Write our graph out to disk:
    print("Writing graph pickle...")
    nx.write_gpickle(graph, f"{out_base_name}.graph.gpickle")
    print("Done!")


##########################################################################################
##########################################################################################
##########################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=f"Processes the refmap files from gffcompare to produce a graph of relationships between gff-reads "
                    f"and the transcripts/genes to which they are aligned.  This graph is then used to produce "
                    f"assignments for each read to a gene and transcript.",
    )

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--gencode_gtf',
                               help='gencode gtf file used in the gffcompare calls for this dataset',
                               required=True)

    requiredNamed.add_argument('--st2-mas-seq-refmap',
                               help='stringtie2 reference <-> MAS-seq query refmap file',
                               required=True)
    requiredNamed.add_argument('--st2-gencode-refmap',
                               help='stringtie2 reference <-> gencode query refmap file',
                               required=True)
    requiredNamed.add_argument('--gencode-mas-seq-refmap',
                               help='Gencode reference <-> MAS-seq query refmap file',
                               required=True)
    requiredNamed.add_argument('--gencode-st2-refmap',
                               help='Gencode reference <-> stringtie2 query refmap file',
                               required=True)

    requiredNamed.add_argument('--st2-mas-seq-tmap',
                               help='stringtie2 reference <-> MAS-seq query tmap file',
                               required=True)
    requiredNamed.add_argument('--st2-gencode-tmap',
                               help='stringtie2 reference <-> gencode query tmap file',
                               required=True)
    requiredNamed.add_argument('--gencode-mas-seq-tmap',
                               help='Gencode reference <-> MAS-seq query tmap file',
                               required=True)
    requiredNamed.add_argument('--gencode-st2-tmap',
                               help='Gencode reference <-> stringtie2 query tmap file',
                               required=True)

    requiredNamed.add_argument('-o', '--out-base-name',
                               help='Base name of the output files that will be produced.',
                               required=True)

    args = parser.parse_args()
    main(args.gencode_gtf,
         args.st2_gencode_refmap, args.st2_mas_seq_refmap, args.gencode_st2_refmap, args.gencode_mas_seq_refmap,
         args.st2_gencode_tmap, args.st2_mas_seq_tmap, args.gencode_st2_tmap, args.gencode_mas_seq_tmap,
         args.out_base_name)

