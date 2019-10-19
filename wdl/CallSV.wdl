version 1.0

import "Structs.wdl"

task PBSV {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fai
        File tandem_repeat_bed

        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GB")) + 10

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        pbsv \
          discover \
          --tandem-repeats \
          ~{tandem_repeat_bed} \
          ~{bam} \
          ~{output_prefix}.svsig.gz

        pbsv \
          call \
          --num-threads ${num_core} \
          ~{ref_fasta} \
          ~{output_prefix}.svsig.gz \
          ~{output_prefix}.pbsv.vcf
    >>>

    output {
        # save both svsig and vcf
        #File raw_signal = "~{output_prefix}.svsig.gz"
        File variants = "~{output_prefix}.pbsv.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4, 
        mem_gb:             26, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/biocontainers/pbsv:2.2.2--0"
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

task Sniffles {
    input {
        File bam
        File bai

        Int? min_read_support = 3

        String sample_name
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GB")) + 2

    command <<<
        set -euxo pipefail

        # The following choices of parameters follow the descriptions in 
        # Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome
        # TODO: we can generalize it later
        sniffles \
          -m ~{bam} \
          -v ~{output_prefix}.sniffles.pre.vcf \
          --min_support ~{min_read_support} \
          --genotype \
          --report_seq \
          --skip_parameter_estimation

        cat ~{output_prefix}.sniffles.pre.vcf | sed 's/FORMAT\t\/cromwell_root.*.bam/FORMAT\t~{sample_name}/' > ~{output_prefix}.sniffles.vcf
    >>>

    output {
        File variants = "~{output_prefix}.sniffles.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4, 
        mem_gb:             15, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/biocontainers/sniffles:1.0.11--hdbcaa40_1"
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

task CompressAndIndex {
    input {
        File vcf

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(vcf)
    Int disk_size = 2*ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        bgzip -c ~{vcf} > ~{basename}.gz
        tabix -p vcf ~{basename}.gz
    >>>

    output {
        File variants = "~{basename}.gz"
        File variants_tbi = "~{basename}.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "kgarimella/lr-align:0.01.18"
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
