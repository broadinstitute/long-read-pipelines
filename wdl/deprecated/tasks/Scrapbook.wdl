version 1.0

task CoverageTrack {
    input {
        File bam
        File bai
        String chr
        String start
        String end

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        samtools depth -a ~{bam} -r ~{chr}:~{start}-~{end} | bgzip > ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
        tabix -p bed ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
    >>>

    output {
        File coverage = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz"
        File coverage_tbi = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task FilterMQ0Reads {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))
    String prefix = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        samtools view -q 1 -b ~{bam} > ~{prefix}.no_mq0.bam
        samtools index ~{prefix}.no_mq0.bam
    >>>

    output {
        File no_mq0_bam = "~{prefix}.no_mq0.bam"
        File no_mq0_bai = "~{prefix}.no_mq0.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task ComputeBedCoverage {
    input {
        File bam
        File bai
        File bed
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(bed, "GB"))

    command <<<
        set -euxo pipefail

        bedtools coverage -b ~{bed} -a ~{bam} -nobuf | gzip > ~{prefix}.txt.gz
        zcat ~{prefix}.txt.gz | awk '{ sum += sprintf("%f", $15*$16) } END { printf("%f\n", sum) }' > ~{prefix}.count.txt
    >>>

    output {
        File coverage = "~{prefix}.txt.gz"
        Float counts = read_float("~{prefix}.count.txt")
        File counts_file = "~{prefix}.count.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task BamToBed {
    input {
        File bam
        File bai

        RuntimeAttr? runtime_attr_override
    }

    String bed = basename(bam, ".bam") + ".bed"
    Int disk_size = 4*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        bedtools bamtobed -i ~{bam} > ~{bed}
    >>>

    output {
        File bedfile = bed
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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
