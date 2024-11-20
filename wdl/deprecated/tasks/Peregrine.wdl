version 1.0

##########################################################################################
# This pipeline assembles sequences an PacBio HiFi BAM using Peregrine
#  https://github.com/cschin/Peregrine
# and aligns the assembly to a reference, followed by variant calling using paftools
#  https://github.com/lh3/minimap2/tree/master/misc#asmvar
##########################################################################################

import "../../structs/Structs.wdl"

workflow Peregrine {
    input {
        File ref_fasta
        File bam
        String sample_name
        String output_prefix
    }

    parameter_meta {
        bam           : "CCS BAM file containing sequences to be assembled"
        ref_fasta     : "reference to which the draft assembly will be aligned to"
        sample_name   : "name of the sample in the BAM"
        output_prefix : "prefix to output VCF"
    }

    call Assemble {
        input:
            bam = bam,
    }

    call AlignAsPAF {
        input:
            ref_fasta = ref_fasta,
            final_fa = Assemble.final_fa
    }

    call CallWithPaftools {
        input:
            ref_fasta = ref_fasta,
            paf = AlignAsPAF.paf,
            sample_name = sample_name,
            output_prefix = output_prefix
    }

    output {
        File final_fa = Assemble.final_fa
        File paf      = AlignAsPAF.paf
        File variants = CallWithPaftools.variants
    }
}

task Assemble {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        df -h .

        samtools fastq ~{bam} | gzip -1 > reads.fastq.gz

        find $PWD -name "*.fastq.gz" | sort > reads.lst

        echo 'yes' | pg_run.py asm reads.lst 24 24 24 24 24 24 24 24 24 \
            --with-consensus --shimmer-r 3 --best_n_ovlp 8 \
            --output ./pasm

        df -h .
        tree -h
    >>>

    output {
        File final_fa = "pasm/p_ctg_cns.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          48,
        mem_gb:             384,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-peregrine:0.1.6"
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

task AlignAsPAF {
    input {
        File ref_fasta
        File final_fa

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(ref_fasta, "GB") + size(final_fa, "GB"))
    Int num_cpus = 4

    command <<<
        set -euxo pipefail

        minimap2 --paf-no-hit -cx asm20 --cs -r 2k -t ~{num_cpus} ~{ref_fasta} ~{final_fa} | gzip -1 > out.paf.gz
    >>>

    output {
        File paf = "out.paf.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          "~{num_cpus}",
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        3,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-peregrine:0.1.6"
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

task CallWithPaftools {
    input {
        File ref_fasta
        File paf
        String sample_name
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(ref_fasta, "GB") + size(paf, "GB"))
    Int num_cpus = 1

    command <<<
        zcat ~{paf} | sort -k6,6 -k8,8n | /opt/test/minimap2-2.17_x64-linux/k8 /opt/test/minimap2-2.17_x64-linux/paftools.js call -f ~{ref_fasta} -s ~{sample_name} - > ~{output_prefix}.peregrine.vcf
    >>>

    output {
        File variants = "~{output_prefix}.peregrine.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          "~{num_cpus}",
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        3,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-peregrine:0.1.6"
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
