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