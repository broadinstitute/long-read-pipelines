#!/usr/bin/env python

import yaml
import argparse
from collections import defaultdict
from enum import Enum
import pysam
import numpy as np
from collections import Counter
import pickle
import operator
from polyleven import levenshtein
from tqdm import tqdm


READ_TYPE = Enum("READ_TYPE", 'CCS '
                              'CLR ')

UMI_LEN = 10
UMI_TAG = "JX" #"ZU"
FINAL_UMI_TAG = "BX"
UMI_CORR_TAG = "UX"
GENE_TAG = "eq"


def get_read_type(read):
    return READ_TYPE.CCS if read.get_tag('rq') != -1 else READ_TYPE.CLR


def get_read_locus(read):
    return (read.get_tag('CB'), read.get_tag(GENE_TAG))


def get_read_seq(read):
    start, end = read.get_tag('SG').split("cDNA:",1)[1].split(",")[0].split("-") 
    return read.query_sequence.upper()[int(start):int(end)+1]


def get_back_aln_score(read):
    return int(read.get_tag('JB').split("/")[0])


class ReadSnapshot:
    def __init__(self, read):
        self.umi = read.get_tag(UMI_TAG)
        self.type = get_read_type(read)
        self.start = read.reference_start
        self.end = read.reference_end
        sequence = get_read_seq(read)
        self.len = len(sequence)
        self.gc = float(sequence.count('C') + sequence.count('G'))/len(sequence)
        self.name = read.qname


def valid_umi(read, config):
    # checks the deviation of the UMI length
    return abs(len(read.get_tag(UMI_TAG)) - UMI_LEN) <= config.max_umi_delta[get_read_type(read).name]


def valid_gene(read):
    # requires either a MAS or ENSG gene tag
    return ("MAS" in read.get_tag("XG")) or ("ENSG" in read.get_tag("XG"))


def valid_tags(read, config):
    # checks for the presence of required tags
    return  read.has_tag('rq') and \
            read.has_tag('CB') and \
            read.has_tag(UMI_TAG) and \
            read.has_tag(GENE_TAG)


def filter_read(read, config):
    # filters the read based on the final UMI length and back alignment score
    return get_back_aln_score(read) >= config.min_back_aln_score and \
           abs(len(read.get_tag(FINAL_UMI_TAG)) - UMI_LEN) <= config.max_umi_delta_filter[get_read_type(read).name]


def extract_read_groups(input_bam_fname, config):  
    locus2reads = defaultdict(list)
    n_filtered_umi = 0
    n_filtered_gene = 0
    n_valid_reads = 0
    with pysam.AlignmentFile(input_bam_fname, "rb") as input_bam:
        for read in tqdm(input_bam):
            if not valid_tags(read, config):
                continue
            if not valid_gene(read):
                n_filtered_gene += 1
                continue
            if not valid_umi(read, config):
                n_filtered_umi += 1
                continue
            n_valid_reads += 1
            locus = get_read_locus(read)
            locus2reads[locus].append(ReadSnapshot(read))
        print("Number of valid reads: ", n_valid_reads)
        print("Number of filtered by gene: ", n_filtered_gene)
        print("Number of filtered by UMI: ", n_filtered_umi)
    return locus2reads


class UMIConfig:
    def __init__(self, **entries):
        self.__dict__.update(entries)

    def __str__(self):
        s = '\n'.join("{}: {}".format(k, v) for k, v in self.__dict__.items())
        return s


# UMI correction: set cover algorithm

def get_conversion_type(source, target):
    return "CCS" if ((source.type == target.type) and (source.type == READ_TYPE.CCS)) else "CLR"


def can_convert(source, target, config):
    conversion_type = get_conversion_type(source, target)
    edit_dist = levenshtein(source.umi, target.umi, config.edit_distance[conversion_type])
    if edit_dist > config.edit_distance[conversion_type]:
        return False
    delta_len = abs(source.len - target.len)    
    delta_gc = abs(source.gc - target.gc)
    op = operator.and_
    if config.op == "OR":
        op = operator.or_
    return op(delta_len <= config.len_diff[conversion_type], 
              delta_gc <= config.gc_diff[conversion_type])


def get_unique_targets(reads):
    return list(set([r for r in reads]))


def build_graph(reads, config):
    targets = get_unique_targets(reads)
    graph = defaultdict(list)
    target2umi_counts = defaultdict(Counter)
    target2umi_seq = {target_id: target.umi for target_id, target in enumerate(targets)}    
    for read_id, read in enumerate(reads):
        for target_id, target in enumerate(targets):
            if can_convert(read, target, config):
                graph[target_id].append(read_id)
                target2umi_counts[target_id][read.umi] += 1
    return graph, target2umi_counts, target2umi_seq


def min_vertex_cover(read_ids, graph, target2umi_counts, target2umi_seq):
    umi_groups = []
    while len(read_ids) > 0:
        # find the largest group, tie breaking: target matching the max support UMI in group
        max_target_id = max(graph.keys(), key=lambda t: (len(graph[t]), max(target2umi_counts[t], key=target2umi_counts[t].get) == target2umi_seq[t]))
        max_size = len(graph[max_target_id])
        if max_size == 1:
            break
        umi_groups.append(graph[max_target_id])
        # remove reads in the largest group from other groups
        for selected_read_id in graph[max_target_id]:
            for target_id in graph:
                if target_id == max_target_id:
                    continue
                for read_id in graph[target_id]:
                    if selected_read_id == read_id:
                        graph[target_id].remove(read_id)
            read_ids.remove(selected_read_id)
        del graph[max_target_id]
    return umi_groups


def process_reads_at_locus(reads, read2umi, config, debug=False):
    if len(reads) < 2:
        return
    # if all the reads have the same umi, skip correction
    if all(read.umi == reads[0].umi for read in reads):
        return
    graph, target2umi_counts, target2umi_seq = build_graph(reads, config)
    read_ids = list(range(len(reads)))
    umi_groups = min_vertex_cover(read_ids, graph, target2umi_counts, target2umi_seq)
    for group in umi_groups:
        umis = [reads[read_id].umi for read_id in group]
        # pick a umi with maximal support and closest len to UMI_LEN
        umi_max = max(umis, key=lambda t: (umis.count(t), -abs(UMI_LEN-len(t))))
        for read_id in group:
            read2umi[reads[read_id].name] = umi_max


def umi_correction(input_bam_fname, output_bam_fname, filter_bam_fname, config):

    # split reads into groups by locus
    locus2reads = extract_read_groups(input_bam_fname, config)

    # correct reads at each locus
    read2umi = {}
    print("Number of loci: ", len(locus2reads))

    with pysam.AlignmentFile(input_bam_fname, "rb") as input_bam:
        with pysam.AlignmentFile(output_bam_fname, "wb", template=input_bam) as output_bam:
            with pysam.AlignmentFile(filter_bam_fname, "wb", template=input_bam) as filter_bam:

                for locus in tqdm(locus2reads):
                    process_reads_at_locus(locus2reads[locus], read2umi, config)

                # output BAM with corrected UMIs
                for read in tqdm(input_bam):
                    if read.qname in read2umi:
                        read.set_tag(FINAL_UMI_TAG, read2umi[read.qname])
                        read.set_tag(UMI_CORR_TAG, 1)
                    else:
                        read.set_tag(FINAL_UMI_TAG, read.get_tag(UMI_TAG))
                        read.set_tag(UMI_CORR_TAG, 0)
                    if filter_read(read, config):
                        filter_bam.write(read)
                    else:
                        output_bam.write(read)


def main():
    parser = argparse.ArgumentParser(description='UMI correction with set cover')
    parser.add_argument('--input_bam', help='Input BAM file', required=True)
    parser.add_argument('--output_bam', help='Output BAM file with corrected UMIs in the BX tag', required=True)
    parser.add_argument('--filtered_bam', help='BAM file with reads filtered based on UMI len/alignment', required=True)
    parser.add_argument('--config', help='YAML configuration file', required=True)
    args = parser.parse_args()
    with open(args.config) as config_file:
        config = yaml.load(config_file, Loader=yaml.FullLoader)
    umi_config = UMIConfig(**config)
    print(umi_config)
    umi_correction(args.input_bam, args.output_bam, args.filtered_bam, umi_config)


if __name__ == '__main__':
    main()




