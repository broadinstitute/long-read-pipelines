version 1.0

import "Structs.wdl"

task RemoveMasSeqTruncatedReads {
    input {
        File bam

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + ceil(2 * size(bam, "GiB"))

    command <<<
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        /python_scripts/remove_mas_seq_trucated_reads.py ~{bam} ~{prefix}.non_truncated
        samtools index -@$np ~{prefix}.non_truncated.bam
    >>>

    output {
        File non_trucated_bam = "~{prefix}.non_truncated.bam"
        File non_trucated_bai = "~{prefix}.non_truncated.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  3,
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

task AdjustUmiSequenceWithAdapterAlignment {
    # TODO: Move this into Longbow - both the WDL and the code itself.
    meta {
        description : "Extracts a new UMI from each given read by aligning the preceding adapter sequences to the read."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }
    input {
        File bam
        String prefix = "out"

        Int umi_length = 10
        String existing_umi_tag = "ZU"
        String new_umi_tag = "JX"

        String? pre_pre_umi_seq

        String? pre_umi_seq
        String? pre_umi_tag

        String? post_umi_seq
        String? post_umi_tag

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 3*ceil(size(bam, "GB"))

    # NOTE: Mutual exclusivity of arguments is coded in the script itself.
    String pre_pre_umi_seq_arg = if defined(pre_pre_umi_seq)  then " --pre-pre-umi-seq " else ""

    String pre_umi_seq_arg = if defined(pre_umi_seq)  then " --pre-umi-seq " else ""
    String pre_umi_tag_arg = if defined(pre_umi_tag)  then " --pre-umi-tag " else ""

    String post_umi_seq_arg = if defined(post_umi_seq)  then " --post-umi-seq " else ""
    String post_umi_tag_arg = if defined(post_umi_tag)  then " --post-umi-tag " else ""

    command {
        set -e

        python3 /lrma/update_umi_positions_2.py \
            -b ~{bam} \
            -s /dev/null \
            --umi-length ~{umi_length} \
            --existing-umi-tag ~{existing_umi_tag} \
            --new-umi-tag ~{new_umi_tag} \
            ~{pre_pre_umi_seq_arg}~{default="" sep=" --pre-pre-umi-seq " pre_pre_umi_seq} \
            ~{pre_umi_seq_arg}~{default="" sep=" --pre-umi-seq " pre_umi_seq} \
            ~{pre_umi_tag_arg}~{default="" sep=" --pre-umi-tag " pre_umi_tag} \
            ~{post_umi_seq_arg}~{default="" sep=" --post-umi-seq " post_umi_seq} \
            ~{post_umi_tag_arg}~{default="" sep=" --post-umi-tag " post_umi_tag} \
            -o ~{prefix}.umi_adjusted.bam | tee ~{prefix}.log
    }

    output {
        File umi_adjusted_bam = "~{prefix}.umi_adjusted.bam"
        File log = "~{prefix}.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.18"
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

task FilterMasSeqReads {
    input {
        File input_bam
        File input_bai

        Int maxReadLength = 15000
        Int maxEndClipping = 1000

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + ceil(2 * size(input_bam, "GiB")) + ceil(size(input_bai, "GiB"))

    command {
        /gatk/gatk PrintReads \
            -I ~{input_bam} \
            -O ~{prefix}.bam \
            --disable-read-filter WellformedReadFilter \
            --read-filter MappedReadFilter \
            --read-filter MappingQualityNotZeroReadFilter \
            --read-filter NotSecondaryAlignmentReadFilter \
            --read-filter NotSupplementaryAlignmentReadFilter \
            --read-filter ReadLengthReadFilter --max-read-length 15000 \
            --read-filter ExcessiveEndClippedReadFilter --max-clipped-bases 1000

        echo "PWD is:"
        pwd

        echo "PWD List:"
        ls -la

        echo "Outfile list:"
        ls -la ~{prefix}.bam*

        date
    }

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "broadinstitute/gatk:4.2.6.1"
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

task RenameSingleCellBamTagsForMasIsoSeqV0 {
    meta {
        description : "Rename the single-cell tags in MAS-ISO-seq v0 data (CB -> Jp; ZU -> Jq ...)."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam

        String prefix = "tags_renamed"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        time /python_scripts/rename_single_cell_bam_tags.py \
            ~{bam} \
            ~{prefix}.bam

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

task RestoreSingleCellBamTagsForMasIsoSeqV0 {
    meta {
        description : "Restore the single-cell tags in MAS-ISO-seq v0 data (Jp -> Cb; Jq -> ZU ...)."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam

        String prefix = "tags_restored"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        time /python_scripts/restore_single_cell_bam_tags.py \
            ~{bam} \
            ~{prefix}.bam

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