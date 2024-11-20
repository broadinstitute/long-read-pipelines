version 1.0

import "../../structs/Structs.wdl"

workflow CallAssemblyVariants {

    meta {
        description: "Call variants from an assembly using paftools.js"
    }

    parameter_meta {
        asm_fasta:        "haploid assembly"
        ref_fasta:        "reference to which assembly should be aligned"
        participant_name: "participant name"
        prefix:           "prefix for output files"
    }

    input {
        File asm_fasta
        File ref_fasta
        String participant_name
        String prefix
    }

    call AlignAsPAF {
        input:
            ref_fasta = ref_fasta,
            asm_fasta = asm_fasta,
            prefix = prefix
    }

    call Paftools {
        input:
            ref_fasta = ref_fasta,
            paf = AlignAsPAF.paf,
            participant_name = participant_name,
            prefix = prefix
    }

    output {
        File paf = AlignAsPAF.paf
        File paftools_vcf = Paftools.variants
    }
}

task AlignAsPAF {
    input {
        File ref_fasta
        File asm_fasta
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(ref_fasta, "GB") + size(asm_fasta, "GB"))
    Int num_cpus = 4

    command <<<
        set -euxo pipefail

        minimap2 --paf-no-hit -cx asm20 --cs -r 2k -t ~{num_cpus} \
            ~{ref_fasta} ~{asm_fasta} | \
            gzip -1 > ~{prefix}.paf.gz
    >>>

    output {
        File paf = "~{prefix}.paf.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             40,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-asm:0.1.13"
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

task Paftools {
    input {
        File ref_fasta
        File paf
        String participant_name
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(ref_fasta, "GB") + size(paf, "GB"))
    Int num_cpus = 1

    command <<<
        zcat ~{paf} | \
            sort -k6,6 -k8,8n | \
            paftools.js call -f ~{ref_fasta} -s ~{participant_name} - \
            > ~{prefix}.paftools.vcf
    >>>

    output {
        File variants = "~{prefix}.paftools.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-asm:0.1.13"
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
