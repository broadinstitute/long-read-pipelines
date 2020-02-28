version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test-2/wdl/tasks/Structs.wdl"

task SelectReadsFromRegion {
    input {
        File bam
        File bai
        Array[String]+ region

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        # select region(s)
        samtools view -hb ~{bam} ~{sep=" " region} | samtools fastq - | gzip -1 > reads.fastq.gz
    >>>

    output {
        File reads = "reads.fastq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.10"
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

task CombineReads {
    input {
        Array[File]+ reads

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        zcat ~{sep=" " reads} | python3 /usr/local/bin/cat_as_fasta.py | gzip -1 > combined.reads.fastq.gz
    >>>

    output {
        File reads = "combined.reads.fastq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.10"
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

task CorrectAndTrimReadsWithCanu {
    input {
        File reads
        String target_size
        String platform
        Boolean? is_corrected
        String prefix

        Int? max_coverage
        RuntimeAttr? runtime_attr_override
    }

    Int readSamplingCoverage = select_first([max_coverage, 500])
    Boolean correct = select_first([is_corrected, false])
    String data_type = "-" + (if platform == "PACBIO" then "pacbio" else "nanopore") + "-" + (if correct then "corrected" else "raw")
    Int disk_size = 4*ceil(size(reads, "GB"))
    Int cpus = 8

    command <<<
        set -euxo pipefail

        # correct and trim reads
        canu -correct -p ~{prefix} -d ./ genomeSize=~{target_size} readSamplingCoverage=~{readSamplingCoverage} ~{data_type} ~{reads}
        canu -trim -p ~{prefix} -d ./ genomeSize=~{target_size} readSamplingCoverage=~{readSamplingCoverage} ~{data_type} ~{prefix}.correctedReads.fasta.gz
    >>>

    output {
        File report          = "~{prefix}.report"
        File corrected_reads = "~{prefix}.correctedReads.fasta.gz"
        File trimmed_reads   = "~{prefix}.trimmedReads.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.10"
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

task AssembleReadsWithCanu {
    input {
        File reads
        String target_size
        String platform
        Boolean? is_corrected
        String prefix

        Int? min_coverage
        Int? max_coverage
        RuntimeAttr? runtime_attr_override
    }

    Int stopOnLowCoverage = select_first([min_coverage, 0])
    Int readSamplingCoverage = select_first([max_coverage, 500])
    Boolean correct = select_first([is_corrected, false])
    String data_type = "-" + (if platform == "PACBIO" then "pacbio" else "nanopore") + "-" + (if correct then "corrected" else "raw")

    String canutask = if correct then "-assemble" else ""
    String errorrate = if correct then "" else "correctedErrorRate=0.105"

    Int disk_size = 4*ceil(size(reads, "GB"))
    Int cpus = 8

    command <<<
        set -euxo pipefail

        # assemble target
        canu ~{canutask} -p ~{prefix} -d ./ genomeSize=~{target_size} stopOnLowCoverage=~{stopOnLowCoverage} readSamplingCoverage=~{readSamplingCoverage} ~{errorrate} ~{data_type} ~{reads}
    >>>

    output {
        File report          = "~{prefix}.report"

        File contigs_fasta   = "~{prefix}.contigs.fasta"
        File unassembled     = "~{prefix}.unassembled.fasta"
        File unitigs_fasta   = "~{prefix}.unitigs.fasta"

        File contigs_layout  = "~{prefix}.contigs.layout"
        File unitigs_layout  = "~{prefix}.unitigs.layout"
        File unitigs_bed     = "~{prefix}.unitigs.bed"

        File contigs_gfa     = "~{prefix}.contigs.gfa"
        File unitigs_gfa     = "~{prefix}.unitigs.gfa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.10"
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

task AlignContigs {
    input {
        File contigs
        File ref_fasta
        String SM
        String ID
        String PL
        Boolean? is_corrected
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Boolean correct = select_first([is_corrected, false])
    String map_arg = if (PL == "ONT") then "map-ont" else "map-pb"
    String correction_arg = if (correct) then "asm20" else map_arg

    Int cpus = 4
    Int disk_size = ceil(size(ref_fasta, "GB")) + 4*ceil(size(contigs, "GB"))

    String aligned_name = "~{prefix}.aligned.bam"

    command <<<
        set -euxo pipefail

        awk '{ print $1 }' ~{contigs} | minimap2 -ayYL --cs --MD --eqx -x ~{correction_arg} -R '@RG\tID:~{ID}\tSM:~{SM}\tPL:~{PL}' -t ~{cpus} ~{ref_fasta} - | samtools sort -@~{cpus} -m4G -o ~{aligned_name} -
        samtools index ~{aligned_name}
    >>>

    output {
        File aligned_bam = "~{aligned_name}"
        File aligned_bai = "~{aligned_name}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.10"
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

task CallHaploidVariants {
    input {
        File bam
        File bai
        File ref_fasta
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 2
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(ref_fasta, "GB"))

    command <<<
        set -euxo pipefail

        bcftools mpileup -Ou -f ~{ref_fasta} ~{bam} | bcftools call --ploidy 1 -mv -Ov -o ~{prefix}.calls.vcf
    >>>

    output {
        File calls = "~{prefix}.calls.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.10"
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
