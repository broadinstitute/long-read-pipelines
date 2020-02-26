version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test/wdl/tasks/Structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test/wdl/tasks/Finalize.wdl" as FF

workflow AssemblyMetrics {
    input {
        File asm
        String asm_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File gff
        File trf

        String? gcs_output_dir
    }

    call ComputeExonAccuracy {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            gff = gff,
            asm = asm,
    }

    call Utils.BamToTable {
        input:
            bam = ComputeExonAccuracy.aligned_exons_bam,
            prefix = asm_name + ".aligned_exons.table.txt"
    }

    output {
        File aligned_exon_table = BamToTable.table
    }
}

task ComputeExonAccuracy {
    input {
        File ref_fasta
        File ref_fasta_fai
        File gff
        File asm

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(ref_fasta, "GiB") + size(ref_fasta_fai, "GiB") + size(gff, "GiB") + size(asm, "GiB"))

    command <<<
        set -euxo pipefail

        awk '{ if ($3 == "exon") print $1 ":" $4 "-" $5 }' ~{gff} > exons.txt
        samtools faidx ~{ref_fasta} -r exons.txt -n 1000000 > exons.fasta

        samtools faidx ~{asm}
        minimap2 -ayYL --MD --eqx -x asm20 ~{asm} exons.fasta | samtools sort | samtools view -b > aligned_exons.bam
        samtools index aligned_exons.bam
    >>>

    output {
        File exons_txt = "exons.txt"
        File exons_fa  = "exons.fasta"
        File aligned_exons_bam = "aligned_exons.bam"
        File aligned_exons_bai = "aligned_exons.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-metrics:0.01.07"
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
