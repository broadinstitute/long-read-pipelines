version 1.0

import "Structs.wdl"

task PhaseCommonVariants {
    input {
        File input_vcf
        Float filter_maf = 0.005
        String interval

        Int num_cpus = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size(input_vcf, "GB"))
    String out_bcf = "common.phased_" + sub(interval, "[:-]", "_") + ".bcf"

    command <<<
        set -euxo pipefail

        SHAPEIT5_phase_common_static \
            --input ~{input_vcf} \
            --filter-maf ~{filter_maf} \
            --output ~{out_bcf} \
            --region ~{interval} \
            --thread ~{num_cpus}

        bcftools index ~{out_bcf} --threads ~{num_cpus}
    >>>

    output {
        File common_phased_shard_bcf = "~{out_bcf}"
        File common_phased_shard_tbi = "~{out_bcf}.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             4*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-shapeit5:1.0.0-beta"
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

task LigatePhasedCommonVariants {
    input {
        Array[File] phased_shard_bcfs
        Array[File] phased_shard_tbis

        String prefix = "out"

        Int num_cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size([phased_shard_bcfs, phased_shard_tbis], "GB"))
    String out_bcf = "~{prefix}.bcf"

    command <<<
        set -euxo pipefail

        SHAPEIT5_ligate_static \
            --input ~{write_lines(phased_shard_bcfs)} \
            --output ~{out_bcf} \
            --thread ~{num_cpus} \
            --index
    >>>

    output {
        File scaffold_bcf = "~{out_bcf}"
        File scaffold_tbi = "~{out_bcf}.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-shapeit5:1.0.0-beta"
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

task PhaseRareVariants {
    input {
        File input_bcf
        File scaffold_bcf
        String input_region
        String scaffold_region

        Int num_cpus = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size([input_bcf, scaffold_bcf], "GB"))
    String out_bcf = "rare.phased_" + sub(input_region, "[:-]", "_") + ".bcf"

    command <<<
        set -euxo pipefail

        SHAPEIT5_phase_rare_static \
            --input-plain ~{input_bcf} \
            --scaffold ~{scaffold_bcf} \
            --output ~{out_bcf} \
            --input-region ~{input_region} \
            --scaffold-region ~{scaffold_region} \
            --thread ~{num_cpus}

        bcftools index ~{out_bcf} --threads ~{num_cpus}
    >>>

    output {
        File rare_phased_shard_bcf = "~{out_bcf}"
        File rare_phased_shard_tbi = "~{out_bcf}.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             4*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-shapeit5:1.0.0-beta"
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

task ConcatenateVariants {
    input {
        Array[File] shard_bcfs
        Array[File] shard_tbis

        String prefix
        Int num_cpus = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size(shard_bcfs, "GB"))
    String out_bcf = "~{prefix}.bcf"

    command <<<
        set -euxo pipefail

        bcftools concat \
            --naive \
            -f ~{write_lines(shard_bcfs)} \
            -o ~{out_bcf} \
            --threads ~{num_cpus}

        bcftools index ~{out_bcf} --threads ~{num_cpus}
    >>>

    output {
        File full_bcf = "~{out_bcf}"
        File full_tbi = "~{out_bcf}.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             4*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-shapeit5:1.0.0-beta"
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
