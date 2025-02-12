version 1.0

import "../../structs/Structs.wdl"
import "../Transcriptomics/MASSeq.wdl" as MAS

struct LongbowModelParams {
    Int umi_length
    String? pre_pre_umi_seq
    String? pre_umi_seq
    String? pre_umi_tag
    String? post_umi_seq
    String? post_umi_tag
}

workflow Process {

    meta {
        description: "Process a BAM file using Longbow"
    }
    parameter_meta {
        bam: "BAM file to process"
        prefix: "Prefix for output files"
        model: "Longbow model to use for processing"
        barcode_allow_list: "File containing a list of barcodes to allow"
        barcode_tag: "Tag containing the barcode"
        corrected_tag: "Tag to use for corrected barcode"
        shard_width: "Width of shards to use for processing"
        same_barcode_per_read: "Whether to assume that all reads in a pair have the same barcode"
    }

    input {
        File bam

        String prefix = "out"
        String? model

        File? barcode_allow_list
        String? barcode_tag
        String? corrected_tag

        Int shard_width = 25

        Boolean same_barcode_per_read = false
    }

    if (!defined(model)) { call Peek as t_01_Peek { input: bam = bam } }
    String lbmodel = select_first([model, t_01_Peek.model])

    # Make a lookup table for our models for now.
    # TODO: Once the UMI adjustment is part of longbow, remove this (it's redundant to Longbow):
    LongbowModelParams umi_params_mas15                = { 'umi_length': 10, 'pre_pre_umi_seq': 'TCTACACGACGCTCTTCCGATCT', "pre_umi_tag": "CB", "post_umi_seq": "TTTCTTATATGGG" }
    LongbowModelParams umi_params_mas15threeP          = { 'umi_length': 12, 'pre_pre_umi_seq': 'TCTACACGACGCTCTTCCGATCT', "pre_umi_tag": "CB", "post_umi_seq": "TTTTTTTTTTTTT" }
    LongbowModelParams umi_params_mas15BulkWithIndices = { 'umi_length': 10, "pre_umi_seq": 'TCTACACGACGCTCTTCCGATCT', "post_umi_seq": "TTTCTTATATGGG" }
    LongbowModelParams umi_params_mas10                = { 'umi_length': 10, 'pre_pre_umi_seq': 'TCTACACGACGCTCTTCCGATCT', "pre_umi_tag": "CB", "post_umi_seq": "TTTCTTATATGGG" }
    LongbowModelParams umi_params_slide_seq            = { 'umi_length': 9, 'pre_pre_umi_seq': 'TCTTCAGCGTTCCCGAGA', "pre_umi_tag": "X1", "post_umi_seq": "TTTTTTTTTTTTT" }
    LongbowModelParams umi_params_mas15teloprimev2     = { 'umi_length': 10, 'pre_pre_umi_seq': 'TCTACACGACGCTCTTCCGATCT', "pre_umi_tag": "CB", "post_umi_seq": "TTTCTTATATGGG" }

    Map[String, LongbowModelParams] longbow_umi_adjustment_params = {
        'mas_15_sc_10x5p_single_none': umi_params_mas15,
        'mas_15_sc_10x3p_single_none': umi_params_mas15threeP,
        'mas_15_bulk_10x5p_single_internal': umi_params_mas15BulkWithIndices,
        'mas_15_spatial_slide-seq_single_none': umi_params_slide_seq,
        'mas_10_sc_10x5p_single_none': umi_params_mas10,
        'mas_15_bulk_teloprimeV2_single_none': umi_params_mas15teloprimev2
    }

    call Annotate as t_02_Annotate { input: prefix = prefix, bam = bam, model = lbmodel }
    call Filter as t_03_Filter { input: prefix = prefix, bam = t_02_Annotate.annotated_bam }
    call Segment as t_04_Segment  { input: prefix = prefix, bam = t_03_Filter.filtered_bam }

    # Here we remove "truncated" reads, which are the final reads of an array that are bounded by the end of the
    # cDNA itself, rather than a known segment delimiter:
    call MAS.RemoveMasSeqTruncatedReads as t_05_RemoveMasSeqTruncatedReads { input: prefix = prefix, bam = t_04_Segment.segmented_bam }

    call PadUMI  as t_07_PadUMI {
        input:
            prefix = prefix + "_umi_padded",
            bam = t_05_RemoveMasSeqTruncatedReads.non_trucated_bam,
            model = lbmodel
    }

    # Only call CBC code if we have a single-cell library:
    if (lbmodel != "mas_15_bulk_10x5p_single_internal" && lbmodel != "mas_15_bulk_10x5p_single_internal") {
        call PadCBC  as t_08_PadCBC {
            input: prefix = prefix + "_umi_padded_cbc_padded",
                bam = t_07_PadUMI.padded_umi_bam,
                model = lbmodel,
                barcode_tag = barcode_tag,
        }
        call Correct as t_09_Correct {
            input: prefix = prefix + "_umi_padded_cbc_padded_corrected",
                bam = t_08_PadCBC.padded_cbc_bam,
                model = lbmodel,
                barcode_allow_list = barcode_allow_list,
                barcode_tag = barcode_tag,
                corrected_tag = corrected_tag
        }
    }

    File bam_for_umi_adjustment = select_first([t_09_Correct.corrected_bam, t_07_PadUMI.padded_umi_bam])

    call MAS.AdjustUmiSequenceWithAdapterAlignment as t_10_AdjustUmiSequenceWithAdapterAlignment {
        input:
            prefix = prefix + "_array_elements_CBC_corrected_UMI_adjusted",
            bam = bam_for_umi_adjustment,
            umi_length = longbow_umi_adjustment_params[lbmodel].umi_length,
            pre_pre_umi_seq = longbow_umi_adjustment_params[lbmodel].pre_pre_umi_seq,
            pre_umi_seq = longbow_umi_adjustment_params[lbmodel].pre_umi_seq,
            pre_umi_tag = longbow_umi_adjustment_params[lbmodel].pre_umi_tag,
            post_umi_seq = longbow_umi_adjustment_params[lbmodel].post_umi_seq,
            post_umi_tag = longbow_umi_adjustment_params[lbmodel].post_umi_tag,
    }

    # Only call CBC code if we have a single-cell library:
    if (lbmodel != "mas_15_bulk_10x5p_single_internal" && lbmodel != "mas_15_bulk_10x5p_single_internal") {
        # Merge our correction stats so we can have a record of them for later:
        call AggregateCorrectLogStats {
            input:
                longbow_correct_log_files = select_all([t_09_Correct.log]),
                out_name = prefix + "_longbow_correct_stats.txt"
        }
    }

    call Extract { input: prefix = prefix, bam = t_10_AdjustUmiSequenceWithAdapterAlignment.umi_adjusted_bam }

    output {
        # Output reads:
        File annotated_bam = t_02_Annotate.annotated_bam
        File segmented_bam = t_04_Segment.segmented_bam

        File filtered_bam = t_03_Filter.filtered_bam
        File filter_failed_bam = t_03_Filter.filter_failed_bam

        File? corrected_bam = t_10_AdjustUmiSequenceWithAdapterAlignment.umi_adjusted_bam
        File? uncorrectable_bam = t_09_Correct.uncorrected_bam

        File extracted_bam = Extract.extracted_bam

        # Output stats / logs:
        File? correct_stats = AggregateCorrectLogStats.stats
        File? correct_log = t_09_Correct.log
        File umi_adjustment_log = t_10_AdjustUmiSequenceWithAdapterAlignment.log
    }
}

task Peek {
    input {
        File bam
        Int n = 100

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow peek -n ~{n} -o model.txt ~{bam}
    >>>

    output {
        String model = read_string("model.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             64,  # TODO: Identify and fix `corrupted double-linked list` issue.  Lots of memory is a bad bandaid.
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task Annotate {
    input {
        File bam
        String model

        String prefix = "out"
        Int num_cpus = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow annotate -m ~{model} -t ~{num_cpus} -o ~{prefix}.annotated.bam ~{bam}
    >>>

    output {
        File annotated_bam = "~{prefix}.annotated.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task Filter {
    input {
        File bam
        File? bam_pbi

        String? model

        Int num_cpus = 2
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    String model_arg = if defined(model) then " --model " else ""
    String pbi_arg = if defined(bam_pbi) then " --pbi " else ""
    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow filter \
            -v INFO \
            ~{pbi_arg}~{default="" bam_pbi} \
            ~{model_arg}~{default="" model} \
            ~{bam} \
            -x ~{prefix}_longbow_filter_failed.bam \
            -o ~{prefix}_longbow_filter_passed.bam 2> >(tee longbow_filter_log.txt >&2) # Get log data from stderr and reprint to stderr
    >>>

    output {
        File filtered_bam = "~{prefix}_longbow_filter_passed.bam"
        File filter_failed_bam = "~{prefix}_longbow_filter_failed.bam"
        File log = "longbow_filter_log.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task Segment {
    input {
        File bam

        Int num_cpus = 2
        String prefix = "out"

        String? model

        RuntimeAttr? runtime_attr_override
    }

    String model_arg = if defined(model) then " --model " else ""
    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow segment -v INFO ~{model_arg}~{default="" model} ~{bam} -o ~{prefix}.segmented.bam
    >>>

    output {
        File segmented_bam = "~{prefix}.segmented.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task Extract {

    input {
        File bam
        Int num_cpus = 2
        Int base_padding = 2

        Int? start_offset         # For mas15: 16+10
        String? leading_adapter   # For mas15: "10x_Adapter"
        String? trailing_adapter  # For mas15: "Poly_A"

        String prefix = "out"

        File? bam_pbi

        RuntimeAttr? runtime_attr_override
    }

    String pbi_arg = if defined(bam_pbi) then " --pbi " else ""
    String start_offset_arg = if defined(start_offset) then " --start-offset " else ""
    String leading_adapter_arg = if defined(leading_adapter) then " --leading-adapter " else ""
    String trailing_adapter_arg = if defined(trailing_adapter) then " --trailing-adapter " else ""

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow extract \
            -v INFO \
            --base-padding ~{base_padding} \
            ~{start_offset_arg}~{default="" start_offset} \
            ~{leading_adapter_arg}~{default=""  leading_adapter} \
            ~{trailing_adapter_arg}~{default=""  trailing_adapter} \
            ~{pbi_arg}~{default="" bam_pbi} \
            ~{bam} \
            -o ~{prefix}.extracted.bam 2> >(tee longbow_extract_log.txt >&2) # Get log data from stderr and reprint to stderr
    >>>

    output {
        File extracted_bam = "~{prefix}.extracted.bam"
        File log = "longbow_extract_log.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task PadUMI
{
    input {
        File bam
        String model

        String prefix = "out"
        Int padding = 2

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(bam, "GB"))

    Map[String, String] umi_tag = {
        "mas_15_sc_10x5p_single_none":          "ZU",
        "mas_15_sc_10x3p_single_none":          "ZU",
        "mas_15_bulk_10x5p_single_internal":    "ZU",
        "mas_10_sc_10x5p_single_none":          "ZU",
        "mas_15_spatial_slide-seq_single_none": "ZU",
        "mas_15_bulk_teloprimeV2_single_none":  "ZU",
    }

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow pad --model ~{model} -v INFO --barcode-tag ~{umi_tag[model]} -e ~{padding} -o tmp.bam -n ~{umi_tag[model]} ~{bam}

        samtools sort tmp.bam -o ~{prefix}.umi_padded_~{padding}.bam
    >>>

    output {
        File padded_umi_bam = "~{prefix}.umi_padded_~{padding}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task PadCBC
{
    input {
        File bam
        String model

        String prefix = "out"
        Int padding = 2

        String? barcode_tag

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(bam, "GB"))

    Map[String, String] barcode_tags = {
        "mas_15_sc_10x5p_single_none":          "CR",
        "mas_15_sc_10x3p_single_none":          "CR",
        "mas_10_sc_10x5p_single_none":          "CR",
        "mas_15_bulk_teloprimeV2_single_none":  "BC",
    }

    String final_barcode_tag = select_first([barcode_tag,barcode_tags[model]])

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow pad --model ~{model} -v INFO --barcode-tag ~{final_barcode_tag} -e ~{padding} -o tmp.bam -n ~{final_barcode_tag} ~{bam}

        samtools sort tmp.bam -o ~{prefix}.cbc_padded_~{padding}.bam
    >>>

    output {
        File padded_cbc_bam = "~{prefix}.cbc_padded_~{padding}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task TagFix
{
    input {
        File bam

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        source /longbow/venv/bin/activate
        longbow tagfix -t${np} -v INFO -o tmp.bam ~{bam}

        samtools sort -@${np} tmp.bam -o ~{prefix}.alignment_tags_fixed.bam
    >>>

    output {
        File tag_fixed_bam = "~{prefix}.alignment_tags_fixed.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task Correct {
    input {
        File bam
        String model

        Int ccs_lev_dist_threshold = 2
        Int clr_lev_dist_threshold = 2

        File? barcode_allow_list
        String? barcode_tag
        String? corrected_tag

        File? barcode_freq_list

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    String barcode_freq_arg = if defined(barcode_freq_list) then " --barcode-freqs " else ""

    Int disk_size = 10 * ceil(size(bam, "GB"))

    Map[String, String] allowlists = {
        "mas_15_sc_10x5p_single_none":          "/longbow/resources/barcodes/cellranger/737K-august-2016.txt",
        "mas_15_sc_10x3p_single_none":          "/longbow/resources/barcodes/cellranger/3M-february-2018.txt.gz",
        "mas_15_bulk_10x5p_single_internal":    "/longbow/resources/barcodes/cellranger/737K-august-2016.txt",
        "mas_10_sc_10x5p_single_none":          "/longbow/resources/barcodes/cellranger/737K-august-2016.txt",
        "mas_15_bulk_teloprimeV2_single_none":  "/longbow/resources/barcodes/indexes/tpv2.indices.txt",
    }

    Map[String, String] machine_memory = {
        "mas_15_sc_10x5p_single_none":          32,
        "mas_15_sc_10x3p_single_none":          64,
        "mas_15_bulk_10x5p_single_internal":    32,
        "mas_10_sc_10x5p_single_none":          32,
        "mas_15_bulk_teloprimeV2_single_none":  32,
    }

    Map[String, String] barcode_tags = {
        "mas_15_sc_10x5p_single_none":          "CR",
        "mas_15_sc_10x3p_single_none":          "CR",
        "mas_15_bulk_10x5p_single_internal":    "CR",
        "mas_10_sc_10x5p_single_none":          "CR",
        "mas_15_bulk_teloprimeV2_single_none":  "BC",
    }

    Map[String, String] corrected_tags = {
        "mas_15_sc_10x5p_single_none":          "CB",
        "mas_15_sc_10x3p_single_none":          "CB",
        "mas_15_bulk_10x5p_single_internal":    "CB",
        "mas_10_sc_10x5p_single_none":          "CB",
        "mas_15_bulk_teloprimeV2_single_none":  "BC",
    }

    # Resolve allow list and barcode tags based on inputs and model:
    String final_barcode_tag = select_first([barcode_tag,barcode_tags[model]])
    String final_corrected_tag = select_first([corrected_tag,corrected_tags[model]])

    command <<<
        set -euxo pipefail

        # For some reason this the allow list specification isn't working in the case of a supplied list.
        # We have to do some cleverness here to make it work:
        user_specified_allow_list="~{barcode_allow_list}"
        if [[ "${#user_specified_allow_list}" -gt 0 ]] ; then
            allow_list=~{barcode_allow_list}
        else
            allow_list=~{allowlists[model]}
        fi

        # NOTE: We can only use 1 thread here because the index will be built independently on each thread,
        #       and the index takes a LOT of memory.
        longbow correct \
            -t 1 \
            --model ~{model} \
            --allow-list ${allow_list} \
            ~{barcode_freq_arg}~{default="" barcode_freq_list} \
            -v INFO \
            --barcode-tag ~{final_barcode_tag} \
            --corrected-tag ~{final_corrected_tag} \
            --max-hifi-dist ~{ccs_lev_dist_threshold} \
            --max-clr-dist ~{clr_lev_dist_threshold} \
            -o ~{prefix}.corrected.bam \
            --barcode-uncorrectable-bam ~{prefix}.uncorrected_barcodes.bam \
            ~{bam} 2>&1 | tee longbow_correct.~{prefix}.log
    >>>

    output {
        File corrected_bam = "~{prefix}.corrected.bam"
        File uncorrected_bam = "~{prefix}.uncorrected_barcodes.bam"
        File log = "longbow_correct.~{prefix}.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             machine_memory[model],  # MUST be this big because of the symspell barcode index.
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task Stats {
    input {
        File bam

        String prefix = "stats"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow stats -o ~{prefix} ~{bam}
    >>>

    output {
        Array[File] pngs = glob("*.png")
        Array[File] svgs = glob("*.svg")
        File summary = "~{prefix}_summary_stats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task Demultiplex {
    input {
        File bam
        String tag

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4 * ceil(size(bam, "GB"))

    String basename = basename(bam, ".bam")

    # TODO: We should do this with 2 passes so that after we've ID'd the models involved, we go back and can refine the demux.  This is mostly for CLR reads.

    command <<<
        set -euxo pipefail

        longbow demultiplex -d ~{tag} -o ~{basename}
    >>>

    output {
        Array[File] demuxed_bams = glob("~{basename}*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

task AggregateCorrectLogStats
{
    input {
        Array[File] longbow_correct_log_files

        String out_name = "longbow_correct_stats.txt"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(longbow_correct_log_files, "GB"))

    # YES, this SHOULD be a proper tool, but right now it isn't.
    command <<<
python << CODE
import os

stats_dict = dict()
line_key = "STATS: "

for stats_file in ["~{sep='","' longbow_correct_log_files}"]:
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
        if tot is None or tot == 0:
            f.write(f"{k}:{' '*k_spacing} {count}\n")
            print(f"WARNING: tot == {tot}")
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
        boot_disk_gb:       25,
        preemptible_tries:  1,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.40"
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

