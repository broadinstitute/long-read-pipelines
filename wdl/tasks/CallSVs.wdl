version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/Structs.wdl"

workflow CallSVs {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File tandem_repeat_bed
    }

    call PBSV {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            prefix = basename(bam, ".bam")
    }

    call Sniffles {
        input:
            bam = bam,
            bai = bai,
            prefix = basename(bam, ".bam")
    }

    call SVIM {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            prefix = basename(bam, ".bam")
    }

    output {
        File pbsv_vcf = PBSV.vcf
        File sniffles_vcf = Sniffles.vcf
        File svim_vcf = SVIM.vcf
    }
}

task PBSV {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File tandem_repeat_bed

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size([bam, bai, ref_fasta, ref_fasta_fai, tandem_repeat_bed], "GB"))

    # purely experiential
    Int memory = if (ceil(size(bam, "GiB")) > 20) then 96 else 64
    Int cpus = ceil( memory / 6 ) # a range of, approximately, [1,6] ratio between mem/cpu allowed from cloud service provider

    command <<<
        set -euxo pipefail

        pbsv discover --tandem-repeats ~{tandem_repeat_bed} ~{bam} ~{prefix}.svsig.gz
        pbsv call --num-threads ~{cpus} ~{ref_fasta} ~{prefix}.svsig.gz ~{prefix}.pbsv.pre.vcf

        cat ~{prefix}.pbsv.pre.vcf | grep -v -e '^chrM' -e '##fileDate' > ~{prefix}.pbsv.vcf
    >>>

    output {
        File vcf = "~{prefix}.pbsv.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.2"
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

        Int min_read_support = 3
        Int min_read_length = 1000
        Int min_mq = 20

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        SM=`samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g'`

        sniffles -t ~{cpus} -m ~{bam} -v ~{prefix}.sniffles.pre.vcf -s ~{min_read_support} -r ~{min_read_length} -q ~{min_mq} --genotype --report_seq --report_read_strands
        sniffles-filter -v ~{prefix}.sniffles.pre.vcf -m ~{min_read_support} -t DEL INS DUP --strand-support 0.001 -l 50 --min-af 0.10 --max-length 400000 -o ~{prefix}.sniffles.filtered.vcf

        cat ~{prefix}.sniffles.filtered.vcf | sed 's/FORMAT\t\/cromwell_root.*.bam/FORMAT\t'"${SM}"'/' | grep -v -e '^chrM' -e '##fileDate' > ~{prefix}.sniffles.vcf
    >>>

    output {
        File vcf = "~{prefix}.sniffles.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             15,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.2"
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

task SVIM {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GB")) + 10

    command <<<
        set -euo pipefail

        SM=`samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g'`

        svim alignment --sample ${SM} --insertion_sequences --read_names ~{prefix}_svim_files ~{bam} ~{ref_fasta}

        grep -v -e '##fileDate' -e '^chrM' ~{prefix}_svim_files/variants.vcf > ~{prefix}.svim.vcf
    >>>

    output {
        File vcf = "~{prefix}.svim.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.2"
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
