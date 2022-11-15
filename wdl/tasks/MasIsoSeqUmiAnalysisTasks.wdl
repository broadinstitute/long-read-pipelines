version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils
import "Longbow.wdl" as LONGBOW
import "Finalize.wdl" as FF

task RestoreSegCoordsAndOriginalNameToReads
{
    input {
        File bam_file
        File? bam_index
        String prefix = "names_and_coords_restored"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam_file, "GB"))

    command <<<
        set -euxo pipefail
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        python3.7 /scripts/00_restore_seg_coords_and_original_name.py ~{bam_file} ~{prefix}

        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mas-iso-seq-umi-analysis:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CopyEqClassInfoToTag
{
    input {
        File bam_file
        File? bam_index
        File eq_class_tsv

        String eq_class_tag = "eq"
        String gene_tag = "XG"

        String prefix = "reads_with_eq_classes"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam_file, "GB")) + 2*ceil(size(eq_class_tsv, "GB"))

    command <<<
        set -euxo pipefail
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        python3.7 /scripts/01_copy_eq_class_info_to_tag.py ~{bam_file} ~{eq_class_tsv} ~{eq_class_tag} ~{gene_tag} ~{prefix}

        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mas-iso-seq-umi-analysis:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ExtractOptimial3pAdapterToUmiTag
{
    input {
        File bam_file
        File? bam_index

        String prefix = "3p_adapters_tagged_as_umis"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam_file, "GB"))

    command <<<
        set -euxo pipefail
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        python3.7 /scripts/02_extract_optimal_3_prime_to_tag.py ~{bam_file} ~{prefix}

        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"

        File ccs_levs_pickle = "~{prefix}.ccs_levs.pickle"
        File clr_levs_pickle = "~{prefix}.clr_levs.pickle"
        File ccs_sws_pickle = "~{prefix}.ccs_sws.pickle"
        File clr_sws_pickle = "~{prefix}.clr_sws.pickle"

        File rejected_bam_no_threep = "~{prefix}.rejected_no_threep.bam"
        File rejected_bam_low_ssw_score = "~{prefix}.rejected_ssw_score_below_35.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mas-iso-seq-umi-analysis:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SplitCcsAndClrReads
{
    input {
        File bam_file
        File? bam_index

        Float min_ccs_rq = 0.0

        String prefix = "reads"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam_file, "GB"))

    command <<<
        set -euxo pipefail
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        set -e

        t_start=$(date +%s.%N)
        echo "Starting CCS / CLR filtration."

        bamtools filter -tag "rq":">=~{min_ccs_rq}" -in ~{bam_file} -out ~{prefix}.ccs.bam &
        bamtools filter  -tag "rq":"<~{min_ccs_rq}" -in ~{bam_file} -out ~{prefix}.clr.bam &

        echo "Waiting for completion."
        wait
        t_end=$(date +%s.%N)
        t_elapsed=$( echo "scale=4;${t_end} - ${t_start}" | bc )

        echo 'Done!'
        echo "Elapsed time: ${t_elapsed}"


        echo "Indexing output files..."
        t_start=$(date +%s.%N)

        samtools index -@$(echo "${np}/2" | bc) ~{prefix}.ccs.bam &
        samtools index -@$(echo "${np}/2" | bc) ~{prefix}.clr.bam &

        echo "Waiting for completion."
        wait
        echo 'Done!'
        echo "Elapsed time: ${t_elapsed}"
    >>>

    output {
        File ccs_bam = "~{prefix}.ccs.bam"
        File ccs_bai = "~{prefix}.ccs.bam.bai"
        File clr_bam = "~{prefix}.clr.bam"
        File clr_bai = "~{prefix}.clr.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mas-iso-seq-umi-analysis:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CreateSimpleCountMatrixForUmiAnalysis
{
    input {
        File bam_file
        File? bam_index

        File eq_class_tsv

        String prefix = "3p_adapters_tagged_as_umis"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam_file, "GB")) + 4*ceil(size(eq_class_tsv, "GB"))

    command <<<
        set -euxo pipefail

        python3.7 /scripts/05_create_simple_count_matrix.py \
            -b ~{bam_file} \
            --tx-eq-class-assignments ~{eq_class_tsv} \
            -o ~{prefix}.tsv
    >>>

    output {
        File simple_counts_tsv = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mas-iso-seq-umi-analysis:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task UmiCoverForThreePrimeAnalysis
{
    meta {
        description : "Uses a set-cover algorithm to correct UMIs, but only for a bespoke 3' adapter analysis.  This is NOT for regular UMI correction."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }
    input {
        File bam
        File? bam_index

        String prefix = "out"
        Int umi_length = 10

        Int max_ccs_edit_dist = 2
        Int max_clr_edit_dist = 3

        Int max_ccs_length_diff = 50
        Int max_clr_length_diff = 150
        Float max_ccs_gc_diff = 0.05
        Float max_clr_gc_diff = 0.15
        Int max_ccs_umi_length_delta = 3
        Int max_clr_umi_length_delta = 4
        Int max_final_ccs_umi_length_delta = 3
        Int max_final_clr_umi_length_delta = 3
        Int min_back_seg_score = 10

        String umi_tag = "JX"
        String gene_tag = "XG"
        String eq_class_tag = "eq"
        String final_umi_tag = "BX"
        String umi_corrected_tag = "UX"
        String back_alignment_score_tag = "JB"

        Boolean pre_extracted = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB")) + 4*ceil(size(bam_index, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        # Generate a config YAML file:
        rm -f correct_umi.yaml

        echo "edit_distance: {\"CCS\": ~{max_ccs_edit_dist}, \"CLR\": ~{max_clr_edit_dist}}" >> correct_umi.yaml
        echo "len_diff: {\"CCS\": ~{max_ccs_length_diff}, \"CLR\": ~{max_clr_length_diff}}" >> correct_umi.yaml
        echo "gc_diff: {\"CCS\": ~{max_ccs_gc_diff}, \"CLR\": ~{max_clr_gc_diff}}" >> correct_umi.yaml
        echo "op: \"AND\"" >> correct_umi.yaml
        echo "max_umi_delta: {\"CCS\": ~{max_ccs_umi_length_delta}, \"CLR\": ~{max_clr_umi_length_delta}}" >> correct_umi.yaml
        echo "max_umi_delta_filter: {\"CCS\": ~{max_final_ccs_umi_length_delta}, \"CLR\": ~{max_final_clr_umi_length_delta}}" >> correct_umi.yaml
        echo "min_back_aln_score: ~{min_back_seg_score}}" >> correct_umi.yaml

python3.7 - --input_bam ~{bam} --output_bam ~{prefix}.corrected_umis.bam --filtered_bam ~{prefix}.failed_umi_correction.bam --config correct_umi.yaml --pre-extracted << CODE
#!/usr/bin/env python3.7

import sys
import time
import os
import pickle
import yaml
import argparse
from functools import reduce
from collections import defaultdict
from enum import Enum
import pysam
import numpy as np
from collections import Counter
import pickle
import operator
from polyleven import levenshtein
from tqdm import tqdm

DEBUG = False


READ_TYPE = Enum("READ_TYPE", 'CCS '
                              'CLR ')

UMI_LEN = 10
UMI_TAG = "~{umi_tag}"
SEG_TAG = "SG"
FINAL_UMI_TAG = "~{final_umi_tag}"
UMI_CORR_TAG = "~{umi_corrected_tag}"
EQ_CLASS_TAG = "~{eq_class_tag}"

READ_QUALITY_TAG = "rq"
CELL_BARCODE_TAG = "CB"
BACK_ALIGNMENT_SCORE_TAG = "~{back_alignment_score_tag}"
GENE_TAG = "~{gene_tag}"

CODING_REGION_NAME = "cDNA"


def get_read_type(read):
    return READ_TYPE.CCS if read.get_tag(READ_QUALITY_TAG) != -1 else READ_TYPE.CLR


def get_read_locus(read):
    return read.get_tag(CELL_BARCODE_TAG), read.get_tag(EQ_CLASS_TAG)


def get_read_seq(read, pre_extracted):
    if pre_extracted:
        return read.query_sequence.upper()
    else:
        start, end = read.get_tag(SEG_TAG).split(f"{CODING_REGION_NAME}:", 1)[1].split(",")[0].split("-")
        if read.is_reverse:
            # RC reads we need to swap coordinates:
            new_start = len(read.query_sequence) - int(end)
            new_end = len(read.query_sequence) - int(start)
            return read.query_sequence.upper()[int(new_start):int(new_end)+1]
        else:
            return read.query_sequence.upper()[int(start):int(end)+1]


def get_back_aln_score(read):
    return 12
    #return int(read.get_tag(BACK_ALIGNMENT_SCORE_TAG).split("/")[0])


class Target:
    def __init__(self, read_snapshot):
        self.umi = read_snapshot.umi
        self.type = read_snapshot.type

    def __eq__(self, other):
        return self.umi == other.umi and self.type == other.type

    def __hash__(self):
        return hash((self.umi, self.type))

class ReadSnapshot:
    def __init__(self, read, pre_extracted):
        self.umi = read.get_tag(UMI_TAG)
        self.type = get_read_type(read)

        #self.start = read.reference_start
        #self.end = read.reference_end
        #sequence = get_read_seq(read, pre_extracted)
        #self.len = len(sequence)
        #self.gc = float(sequence.count('C') + sequence.count('G'))/len(sequence)

        self.name = read.qname


def valid_umi(read, config):
    # checks the deviation of the UMI length
    return abs(len(read.get_tag(UMI_TAG)) - UMI_LEN) <= config.max_umi_delta[get_read_type(read).name]


def valid_gene(read):
    # requires either a MAS or ENSG gene tag
    return ("MAS" in read.get_tag(GENE_TAG)) or ("ENSG" in read.get_tag(GENE_TAG))


def valid_tags(read, config):
    # checks for the presence of required tags
    return  read.has_tag(READ_QUALITY_TAG) and \
            read.has_tag(CELL_BARCODE_TAG) and \
            read.has_tag(UMI_TAG) and \
            read.has_tag(EQ_CLASS_TAG)


def read_passes_filters(read, config):
    # filters the read based on the final UMI length and back alignment score
    return get_back_aln_score(read) >= config.min_back_aln_score and \
           abs(len(read.get_tag(FINAL_UMI_TAG)) - UMI_LEN) <= config.max_umi_delta_filter[get_read_type(read).name]


def extract_read_groups(input_bam_fname, config, pre_extracted):
    locus2reads = defaultdict(list)
    n_filtered_umi = 0
    n_filtered_gene = 0
    n_valid_reads = 0
    n_total_reads = 0
    with pysam.AlignmentFile(input_bam_fname, "rb", check_sq=False, require_index=False) as input_bam:

        total_reads = None
        if input_bam.has_index():
            idx_stats = input_bam.get_index_statistics()
            unaligned_reads = input_bam.nocoordinate
            aligned_reads = reduce(lambda a, b: a + b, [x.total for x in idx_stats]) if len(idx_stats) > 0 else 0
            total_reads = unaligned_reads + aligned_reads

        for read in tqdm(input_bam, desc="Extracting Read Groups", total=total_reads, unit=" read"):
            n_total_reads += 1
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
            locus2reads[locus].append(ReadSnapshot(read, pre_extracted))
        print("Number of reads: ", n_total_reads)
        print("Number of valid reads: ", n_valid_reads)
        print("Number of filtered by gene: ", n_filtered_gene)
        print("Number of filtered by UMI: ", n_filtered_umi)
    return locus2reads, n_total_reads


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

    # For the use case of creating this figure, we don't need to do these checks below, so we can return here:
    return edit_dist <= config.edit_distance[conversion_type]

    #
    # delta_len = abs(source.len - target.len)
    # delta_gc = abs(source.gc - target.gc)
    # op = operator.and_
    # if config.op == "OR":
    #     op = operator.or_
    # return op(delta_len <= config.len_diff[conversion_type],
    #           delta_gc <= config.gc_diff[conversion_type])


def get_unique_targets(reads):
    # TODO: Here extract targets from reads
    return list(set([Target(r) for r in reads]))


def build_graph(reads, config):
    targets = get_unique_targets(reads)

    #if len(reads) * len(targets) > 1e6:
    #    print(f"BIG Group: {len(reads) * len(targets)} : len(reads): {len(reads)} | len(targets): {len(targets)}")
    #    print("READS:")
    #    for r in reads:
    #        print(f"{r.umi}\t{r.type}\t{r.name}")
    #    print("TARGETS:")
    #    for t in targets:
    #        print(f"{t.umi}\t{t.type}")

    graph = defaultdict(list)
    target2umi_counts = defaultdict(Counter)
    target2umi_seq = {target_id: target.umi for target_id, target in enumerate(targets)}
    if len(reads) * len(targets) < 1e6:
        for read_id, read in enumerate(reads):
            for target_id, target in enumerate(targets):
                if can_convert(read, target, config):
                    graph[target_id].append(read_id)
                    target2umi_counts[target_id][read.umi] += 1
    else:
        for read_id, read in tqdm(enumerate(reads), desc="Building graph", position=1):
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
        if DEBUG:
            print("len(reads) < 2")
        return
    # if all the reads have the same umi, skip correction
    if all(read.umi == reads[0].umi for read in reads):
        if DEBUG:
            print("all(read.umi == reads[0].umi for read in reads)")
        return
    graph, target2umi_counts, target2umi_seq = build_graph(reads, config)
    read_ids = list(range(len(reads)))
    umi_groups = min_vertex_cover(read_ids, graph, target2umi_counts, target2umi_seq)
    for group in umi_groups:
        umis = [reads[read_id].umi for read_id in group]
        # pick a umi with maximal support and closest len to UMI_LEN
        umi_max = max(umis, key=lambda t: (umis.count(t), -abs(UMI_LEN-len(t))))
        for read_id in group:
            # TODO: Break up structure so target is just type + umi and leave read the same.
            read2umi[reads[read_id].name] = umi_max


def umi_correction(input_bam_fname, output_bam_fname, filter_bam_fname, config, pre_extracted):

    # split reads into groups by locus
    if os.path.exists(input_bam_fname + ".locus2reads.pickle"):
        print(f"Loading locus2reads cache: {input_bam_fname}.locus2reads.pickle")
        st = time.time()
        locus2reads = pickle.load(open(input_bam_fname + ".locus2reads.pickle", 'rb'))
        try:
            total_num_reads = pickle.load(open(input_bam_fname + ".total_num_reads.pickle", 'rb'))
        except FileNotFoundError:
            total_num_reads = None
        et = time.time()
        print(f"done. (elapsed: {et - st}s)")
    else:
        locus2reads, total_num_reads = extract_read_groups(input_bam_fname, config, pre_extracted)
        print(f"Writing locus2reads cache: {input_bam_fname}.locus2reads.pickle")
        st = time.time()
        pickle.dump(locus2reads, open(input_bam_fname + ".locus2reads.pickle", 'wb'))
        pickle.dump(total_num_reads, open(input_bam_fname + ".total_num_reads.pickle", 'wb'))
        et = time.time()
        print(f"done. (elapsed: {et - st}s)")

    # correct reads at each locus
    read2umi = {}
    print("Number of loci: ", len(locus2reads))

    with pysam.AlignmentFile(input_bam_fname, "rb", check_sq=False, require_index=False) as input_bam:
        with pysam.AlignmentFile(output_bam_fname, "wb", template=input_bam) as correct_umi_bam:
            with pysam.AlignmentFile(filter_bam_fname, "wb", template=input_bam) as filtered_out_umi_bam:

                for locus in tqdm(locus2reads, desc="Processing each locus", position=0):
                    process_reads_at_locus(locus2reads[locus], read2umi, config)

                # output BAM with corrected UMIs
                for read in tqdm(input_bam, desc="Writing out UMI-corrected reads", unit=" read", total=total_num_reads, position=0):
                    if read.qname in read2umi:
                        read.set_tag(FINAL_UMI_TAG, read2umi[read.qname])
                        read.set_tag(UMI_CORR_TAG, 1)
                    else:
                        read.set_tag(FINAL_UMI_TAG, read.get_tag(UMI_TAG))
                        read.set_tag(UMI_CORR_TAG, 0)

                    if read_passes_filters(read, config):
                        correct_umi_bam.write(read)
                    else:
                        filtered_out_umi_bam.write(read)


def main():
    parser = argparse.ArgumentParser(description='UMI correction with set cover')
    parser.add_argument('--input_bam', help='Input BAM file', required=True)
    parser.add_argument('--output_bam', help=f'Output BAM file with corrected UMIs in the "{FINAL_UMI_TAG}" tag', required=True)
    parser.add_argument('--filtered_bam', help='BAM file with reads filtered based on UMI len/alignment', required=True)
    parser.add_argument('--config', help='YAML configuration file', required=True)
    parser.add_argument('--pre-extracted', help='Whether the input file has been processed with longbow extract (default: %(default))', required=False, default=False, action='store_true')
    args = parser.parse_args()
    with open(args.config) as config_file:
        config = yaml.load(config_file, Loader=yaml.FullLoader)
    umi_config = UMIConfig(**config)
    print(umi_config)
    umi_correction(args.input_bam, args.output_bam, args.filtered_bam, umi_config, args.pre_extracted)

if __name__ == '__main__':
    main()

CODE

        samtools sort -@${np} ~{prefix}.corrected_umis.bam > tmp.bam
        mv tmp.bam ~{prefix}.corrected_umis.bam
        samtools index -@${np} ~{prefix}.corrected_umis.bam
    >>>

    output {
        File umi_corrected_bam = "~{prefix}.corrected_umis.bam"
        File umi_corrected_bam_index = "~{prefix}.corrected_umis.bam.bai"
        File failed_umi_correction_bam = "~{prefix}.failed_umi_correction.bam"
        File cached_read_loci = "~{bam}.locus2reads.pickle"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mas-iso-seq-umi-analysis:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
