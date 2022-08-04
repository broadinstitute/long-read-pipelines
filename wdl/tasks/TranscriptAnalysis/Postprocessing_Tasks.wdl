version 1.0

import "../Structs.wdl"

task CreateCountMatrixFromAnnotatedBam {

    meta {
        description : "Creates a count matrix TSV file from the given annotated bam file.  Bam file must contain tags that indicate the gene/transcript (XG), cell barcode (CB), and umi (BX) of the read."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File annotated_transcriptome_bam

        File? tx_equivalence_class_assignments

        String prefix = "umi_tools_group"

        String umi_tag = "ZU"

        RuntimeAttr? runtime_attr_override
    }

    String tx_eq_class_assignments_arg = if defined(tx_equivalence_class_assignments) then " --tx-eq-class-assignments " else ""

    Int disk_size_gb = 20 + 11*ceil(size(annotated_transcriptome_bam, "GB"))
                               + 2*ceil(size(tx_equivalence_class_assignments, "GB"))

    command <<<
        set -euxo pipefail
        /python_scripts/create_count_matrix_from_annotated_bam.py \
            -b ~{annotated_transcriptome_bam} \
            ~{tx_eq_class_assignments_arg} ~{default="" sep=" --tx-eq-class-assignments " tx_equivalence_class_assignments} \
            --umi-tag ~{umi_tag} \
            -o ~{prefix}.tsv
    >>>

    output {
        File count_matrix = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task AggregateUmiAdjustmentStats
{

    # TODO: FINISHME
    input {
        Array[File] longbow_umi_adjustment_log_files

        String out_name = "longbow_umi_adjustment_stats.txt"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(longbow_umi_adjustment_log_files, "GB"))

    # YES, this SHOULD be a proper tool, but right now it isn't.
    command <<<

        f=~{write_lines(longbow_umi_adjustment_log_files)}

        mv $f THE_UMI_FILE_LIST.txt

python << CODE
import os

stats_dict = dict()
line_key = "STATS: "

with open("THE_UMI_FILE_LIST.txt", 'r') as umi_file_list_file:
    for line in umi_file_list_file:
        stats_file = line.strip()
        with open(stats_file, 'r') as f:
            for line in f:
                if line_key in line:
                    line = line.strip()
                    s = line[line.find(line_key) + len(line_key):]
                    key, remainder = [t.strip() for t in s.split(":")]
                    if "/" in remainder:
                        count = int(remainder[:remainder.find("/")])
                        tot = int(remainder[remainder.find("/")+1:remainder.find(" ")])
                    else:
                        count = int(remainder)
                        tot = None

                    try:
                        c, t = stats_dict[key]
                        if tot is not None:
                            tot += t
                        stats_dict[key] = (count + c, tot)
                    except KeyError:
                        stats_dict[key] = (count, tot)

k_len = 0
for k in stats_dict.keys():
    if len(k) > k_len:
        k_len = len(k)

k_prefix = list(stats_dict.keys())[0]
k_prefix = k_prefix[:k_prefix.find(" ")]
with open("~{out_name}", 'w') as f:
    for k, v in stats_dict.items():

        if not k.startswith(k_prefix):
            f.write("\n")
            k_prefix = k[:k.find(" ")]

        k_spacing = k_len - len(k)

        count, tot = v
        if tot is None:
            f.write(f"{k}:{' '*k_spacing} {count}\n")
        else:
            f.write(f"{k}:{' '*k_spacing} {count}/{tot} ({100.0*count/tot:2.4f}%)\n")

CODE
    >>>

    output {
        File stats = out_name
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.27"
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

task MergeBarcodeCounts {
    meta {
        description : "Merge all counts for each unique barcode in the given TSV file.  Assumes file is unheadered and have two columns: BARCODE COUNT.  Merging performed by adding all COUNTs for each BARCODE."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File barcode_count_tsv
        String prefix = "merged_counts"

        RuntimeAttr? runtime_attr_override
    }

    # 20 gb - baseline storage for safety
    # 1x for the file itself
    # 2x for the results and some wiggle room.
    Int disk_size_gb = 20 + (3 * ceil(size(barcode_count_tsv, "GB")))

    command {
        /python_scripts/merge_barcode_counts.py ~{barcode_count_tsv}
        if [[ "~{prefix}.tsv" != "merged_counts.tsv" ]] ; then
            mv merged_counts.tsv "~{prefix}.tsv"
        fi
    }

    output {
        File merged_counts = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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


task CreateCountMatrixAnndataFromTsv {

    meta {
        description : "Creates a python anndata object from the given countmatrix tsv.  Expects the input to have been generated by CreateCountMatrixFromAnnotatedBam.  The resulting anndata object can be directly read into scanpy for single-cell analysis."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File count_matrix_tsv
        File genome_annotation_gtf_file

        Boolean force_anndata_gencode_overwrite = false

        String prefix = "umi_tools_group"

        File? equivalence_class_definitions

        File? overlap_intervals
        String? overlap_interval_label
        File? gencode_reference_gtf_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 20 + 4*ceil(size(count_matrix_tsv, "GB")) + 4*ceil(size(genome_annotation_gtf_file, "GB")) + 2*ceil(size(equivalence_class_definitions, "GB"))

    String overlap_intervals_arg = if defined(overlap_intervals)  then " --overlap-intervals " else ""
    String overlap_interval_label_arg = if defined(overlap_interval_label) then " --overlap-interval-label " else ""
    String gencode_reference_gtf_file_arg = if defined(gencode_reference_gtf_file) then " --gencode-reference-gtf " else ""

    String force_gencode_overwrite_flag = if force_anndata_gencode_overwrite then " --force-overwrite-gencode-overlaps " else ""

    String eq_class_arg = if defined(equivalence_class_definitions)  then " --eq-class-defs-tsv " else ""

    command <<<
        set -euxo pipefail
        /python_scripts/create_count_matrix_anndata_from_tsv.py \
            -t ~{count_matrix_tsv} \
            -g ~{genome_annotation_gtf_file} \
            ~{overlap_intervals_arg}~{default="" sep=" --overlap-intervals " overlap_intervals} \
            ~{overlap_interval_label_arg}~{default="" sep=" --overlap-interval-label " overlap_interval_label} \
            ~{gencode_reference_gtf_file_arg}~{default="" sep=" --gencode-reference-gtf " gencode_reference_gtf_file} \
            ~{eq_class_arg} ~{default="" sep=" --eq-class-defs-tsv " equivalence_class_definitions} \
            ~{force_gencode_overwrite_flag} \
            -o ~{prefix}
    >>>

    output {
        File transcript_gene_count_anndata_h5ad = "~{prefix}_tx_gene_counts_adata.h5ad"
        Array[File] pickles = glob("*.pickle")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task CreateCountMatrixAnndataFromEquivalenceClasses {

    meta {
        description : "Creates a python anndata object from the given countmatrix tsv and equivalence classes.  Expects the input to have been generated by CreateCountMatrixFromAnnotatedBam.  The resulting anndata object can be directly read into scanpy for single-cell analysis."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File count_matrix_tsv
        File genome_annotation_gtf_file

        File tx_equivalence_class_definitions
        File tx_equivalence_class_assignments
        File gene_equivalence_class_definitions
        File gene_equivalence_class_assignments

        Boolean force_anndata_gencode_overwrite = false

        String prefix = "umi_tools_group"

        File? overlap_intervals
        String? overlap_interval_label
        File? gencode_reference_gtf_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 20 + 4*ceil(size(count_matrix_tsv, "GB"))
                          + 4*ceil(size(genome_annotation_gtf_file, "GB"))
                          + 2*ceil(size(tx_equivalence_class_definitions, "GB"))
                          + 2*ceil(size(tx_equivalence_class_assignments, "GB"))
                          + 2*ceil(size(gene_equivalence_class_definitions, "GB"))
                          + 2*ceil(size(gene_equivalence_class_assignments, "GB"))
                          + 2*ceil(size(gencode_reference_gtf_file, "GB"))
                          + 2*ceil(size(overlap_intervals, "GB"))

    String overlap_intervals_arg = if defined(overlap_intervals)  then " --overlap-intervals " else ""
    String overlap_interval_label_arg = if defined(overlap_interval_label) then " --overlap-interval-label " else ""
    String gencode_reference_gtf_file_arg = if defined(gencode_reference_gtf_file) then " --gencode-reference-gtf " else ""

    String force_gencode_overwrite_flag = if force_anndata_gencode_overwrite then " --force-overwrite-gencode-overlaps " else ""

    command <<<
        set -euxo pipefail
        /python_scripts/create_count_matrix_anndata_from_equivalence_classes.py \
            -t ~{count_matrix_tsv} \
            -g ~{genome_annotation_gtf_file} \
            --tx-eq-class-definitions ~{tx_equivalence_class_definitions} \
            --tx-eq-class-assignments ~{tx_equivalence_class_assignments} \
            --gene-eq-class-definitions ~{gene_equivalence_class_definitions} \
            --gene-eq-class-assignments ~{gene_equivalence_class_definitions} \
            ~{overlap_intervals_arg}~{default="" sep=" --overlap-intervals " overlap_intervals} \
            ~{overlap_interval_label_arg}~{default="" sep=" --overlap-interval-label " overlap_interval_label} \
            ~{gencode_reference_gtf_file_arg}~{default="" sep=" --gencode-reference-gtf " gencode_reference_gtf_file} \
            ~{force_gencode_overwrite_flag} \
            -o ~{prefix}
    >>>

    output {
        File transcript_gene_count_anndata_h5ad = "~{prefix}_tx_gene_counts_adata.h5ad"
        Array[File] pickles = glob("*.pickle")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task SubsetCountsMatrixByGenes {

    meta {
        description : "Subsets a count matrix TSV file to contain only the transcripts from the given list of genes.  Assumes the count matrix was produced by comparison with Gencode (due to data formatting) and that the table is a TSV with samples as rows and transcripts as columns."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File count_matrix_tsv
        Array[String] gene_names
    }

    parameter_meta {
        count_matrix_tsv : "TSV file containing the counts of each transcript expressed in a sample.  (Transcripts in columns.  Samples in rows.  One header row.)"
        gene_names : "Array of gene names for which to keep data from the given count matrix TSV."
    }

    # We're subsetting the file, so we should be able to get away with very little space here.
    # 1x for the file itself
    # 2x for the results and some wiggle room.
    Int disk_size = 3 * ceil(size(count_matrix_tsv, "GB"))

    command {
        /python_scripts/subset_count_matrix_by_gene.py ~{count_matrix_tsv} ~{sep=' ' gene_names}
    }

    output {
        File subset_count_matrix_tsv = "count_matrix_subset_by_gene.tsv"
        File subset_count_matrix_h5ad = "count_matrix_subset_by_gene.h5ad"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 8
    }
}


task QuantifyGffComparison {
    meta {
        description : "Create equivalence classes and gene assignments from a set of gffcompare results."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File genome_gtf

        File st2_gencode_refmap
        File st2_gencode_tmap
        File st2_read_refmap
        File st2_read_tmap
        File gencode_st2_refmap
        File gencode_st2_tmap
        File gencode_read_refmap
        File gencode_read_tmap

        String prefix = "reads_comparison"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        genome_gtf : "Genome annotation GTF file (usually gencode)."
        st2_gencode_refmap : "Refmap file (produced by gffcompare) comparing the stringtie2 discovered transcriptome to the genome reference gtf."
        st2_gencode_tmap : "Tmap file (produced by gffcompare) comparing the stringtie2 discovered transcriptome to the genome reference gtf."
        st2_read_refmap : "Refmap file (produced by gffcompare) comparing the stringtie2 discovered transcriptome to input reads (in GFF format)."
        st2_read_tmap : "Tmap file (produced by gffcompare) comparing the stringtie2 discovered transcriptome to input reads (in GFF format)."
        gencode_st2_refmap : "Refmap file (produced by gffcompare) comparing the genome reference gtf to the stringtie2 discovered transcriptome."
        gencode_st2_tmap : "Tmap file (produced by gffcompare) comparing the genome reference gtf to the stringtie2 discovered transcriptome."
        gencode_read_refmap : "Refmap file (produced by gffcompare) comparing the genome reference gtf to input reads (in GFF format)."
        gencode_read_tmap : "Tmap file (produced by gffcompare) comparing the genome reference gtf to input reads (in GFF format)."
        prefix : "Prefix for ouput file."
    }

    Int disk_size_gb = 10 + 2*ceil(size(genome_gtf, "GB"))
                          + 2*ceil(size(st2_gencode_refmap, "GB"))
                          + 2*ceil(size(st2_gencode_tmap, "GB"))
                          + 2*ceil(size(st2_read_refmap, "GB"))
                          + 2*ceil(size(st2_read_tmap, "GB"))
                          + 2*ceil(size(gencode_st2_refmap, "GB"))
                          + 2*ceil(size(gencode_st2_tmap, "GB"))
                          + 2*ceil(size(gencode_read_refmap, "GB"))
                          + 2*ceil(size(gencode_read_tmap, "GB"))

    command <<<
        time /python_scripts/quantify_gff_reads.py \
            --gencode_gtf ~{genome_gtf} \
            --st2-gencode-refmap ~{st2_gencode_refmap} \
            --st2-mas-seq-refmap ~{st2_read_refmap} \
            --gencode-st2-refmap ~{gencode_st2_refmap} \
            --gencode-mas-seq-refmap ~{gencode_read_refmap} \
            --st2-gencode-tmap ~{st2_gencode_tmap} \
            --st2-mas-seq-tmap ~{st2_read_tmap} \
            --gencode-st2-tmap ~{gencode_st2_tmap} \
            --gencode-mas-seq-tmap ~{gencode_read_tmap} \
            -o ~{prefix}

        # Here so the task changes and we don't get cached:
        echo ""
        echo ""
    >>>

    output {
        File gene_eq_class_labels_file = prefix + ".gene_equivalence_class_lookup.tsv"
        File gene_assignments_file = prefix + ".gene_name_assignments.tsv"
        File tx_equivalence_class_labels_file = prefix + ".equivalence_class_lookup.tsv"
        File tx_equivalence_class_file = prefix + ".equivalence_classes.tsv"
        File graph_gpickle = prefix + ".graph.gpickle"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task CombineEqClassFiles {
    meta {
        description : "Combine equivalence classes and gene assignments from disjoint sets of reads produced by QuantifyGffComparison."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        Array[File] gene_eq_class_definitions
        Array[File] gene_assignment_files

        Array[File] equivalence_class_definitions
        Array[File] equivalence_classes

        String prefix = "combined"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gene_eq_class_definitions : "TSV files containing equivalence class definitions for genes as produced by QuantifyGffComparison.gene_eq_class_labels_file."
        gene_assignment_files : "TSV files containing read -> gene equivalence class assignments as produced by QuantifyGffComparison.gene_assignments_file."
        equivalence_class_definitions : "TSV files containing transcript equivalence class definitions as produced by QuantifyGffComparison.tx_equivalence_class_labels_file."
        equivalence_classes : "TSV files containing read -> transcript equivalence class assignments as produced by QuantifyGffComparison.tx_equivalence_class_file."
        prefix : "Prefix for ouput file."
    }

    Int disk_size_gb = 10 + 2*ceil(size(equivalence_class_definitions, "GB"))
                          + 2*ceil(size(equivalence_classes, "GB"))
                          + 2*ceil(size(gene_eq_class_definitions, "GB"))
                          + 2*ceil(size(gene_assignment_files, "GB"))

    command <<<
        time /python_scripts/combine_tx_equivalence_classes.py \
            ~{sep=" " gene_eq_class_definitions} \
            ~{sep=" " gene_assignment_files} \
            ~{sep=" " equivalence_class_definitions} \
            ~{sep=" " equivalence_classes}

        mv gene_equivalence_class_lookup.tsv ~{prefix}.gene_equivalence_class_lookup.tsv
        mv gene_name_assignments.tsv ~{prefix}.gene_name_assignments.tsv
        mv equivalence_class_lookup.tsv ~{prefix}.equivalence_class_lookup.tsv
        mv equivalence_classes.tsv ~{prefix}.equivalence_classes.tsv
    >>>

    output {
        File combined_gene_eq_class_defs = "~{prefix}.gene_equivalence_class_lookup.tsv"
        File combined_gene_eq_class_assignments = "~{prefix}.gene_name_assignments.tsv"
        File combined_tx_eq_class_defs = "~{prefix}.equivalence_class_lookup.tsv"
        File combined_tx_eq_class_assignments = "~{prefix}.equivalence_classes.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task CopyEqClassInfoToTag {
    meta {
        description : "Copy the gene assignment for each given read into the given tag for each read."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam
        File eq_class_file

        String gene_tag = "XG"
        String eq_class_tag = "eq"

        String prefix = "combined"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2*ceil(size(bam, "GB"))
                          + 4*ceil(size(eq_class_file, "GB"))

    # TODO: Extract this code into its own module / file.
    command <<<
        set -euxo pipefail

python << CODE
import pysam
from tqdm import tqdm

with open(f"~{eq_class_file}", 'r') as f:
    read_gene_dict = dict()
    for line in tqdm(f):
        if line.startswith("#"):
            continue
        read_name, tx_eq_class, gene_assignment = line.strip().split("\t")
        read_gene_dict[read_name] = (tx_eq_class, gene_assignment)

with pysam.AlignmentFile(f"~{bam}", "rb", check_sq=False, require_index=False) as bam_file:
    with pysam.AlignmentFile(f"~{prefix}.bam", "wb", header=bam_file.header) as out_bam_file:
        for read in tqdm(bam_file):
            read.set_tag(f"~{eq_class_tag}", read_gene_dict[read.query_name][0])
            read.set_tag(f"~{gene_tag}", read_gene_dict[read.query_name][1])
            out_bam_file.write(read)
CODE

        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam_out = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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
