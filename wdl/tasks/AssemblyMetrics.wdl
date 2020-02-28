version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils
import "Finalize.wdl" as FF

workflow AssemblyMetrics {
    input {
        Array[File] asm
        Array[String] asm_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File gff
        File trf

        String? gcs_output_dir
    }

    scatter (i in range(length(asm))) {
        call ComputeGeneAccuracy {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                gff = gff,
                asm = asm[i],
        }

        call Utils.BamToTable as GeneTable {
            input:
                bam = ComputeGeneAccuracy.aligned_genes_bam,
                prefix = asm_name[i] + ".aligned_genes.table"
        }

        call ComputeVarExonAccuracy {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                gff = gff,
                asm = asm[i],
        }

        call Utils.BamToTable as ExonTable {
            input:
                bam = ComputeVarExonAccuracy.aligned_var_exons_bam,
                prefix = asm_name[i] + ".aligned_var_exons.table"
        }

        call PlotAsmVsRef {
            input:
                asm = asm[i],
                asm_name = asm_name[i],
                ref_fasta = ref_fasta
        }
    }

    call Quast {
        input:
            asm = asm,
            asm_name = asm_name,
            ref_fasta = ref_fasta,
            gff = gff
    }

    output {
        Array[File] aligned_gene_table = GeneTable.table
    }
}

task ComputeGeneAccuracy {
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

        awk '{ if ($3 == "gene") print $1 ":" $4 "-" $5 }' ~{gff} > genes.txt
        samtools faidx ~{ref_fasta} -r genes.txt -n 1000000 > genes.fasta

        samtools faidx ~{asm}
        minimap2 -ayYL --MD --eqx -x asm20 ~{asm} genes.fasta | samtools sort | samtools view -b > aligned_genes.bam
        samtools index aligned_genes.bam
    >>>

    output {
        File genes_txt = "genes.txt"
        File genes_fa  = "genes.fasta"
        File aligned_genes_bam = "aligned_genes.bam"
        File aligned_genes_bai = "aligned_genes.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm-eval:0.01.02"
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

task ComputeVarExonAccuracy {
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

        grep -v '^#' ~{gff} | awk '{ if ($3 == "gene") print $0 }' | grep -i pfemp1 | grep -v -e pseudogene -e exon | sed 's/;.*//' | sed 's/.*ID=//g' | grep -f /dev/stdin ~{gff} | awk '{ if ($3 == "exon") print $0 }' | awk '{ print $1 ":" $4 "-" $5 }' > var_exons.txt
        samtools faidx ~{ref_fasta} -r var_exons.txt -n 1000000 > var_exons.fasta

        samtools faidx ~{asm}
        minimap2 -ayYL --MD --eqx -x asm20 ~{asm} var_exons.fasta | samtools sort | samtools view -b > aligned_var_exons.bam
        samtools index aligned_var_exons.bam
    >>>

    output {
        File var_exons_txt = "var_exons.txt"
        File var_exons_fa  = "var_exons.fasta"
        File aligned_var_exons_bam = "aligned_var_exons.bam"
        File aligned_var_exons_bai = "aligned_var_exons.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm-eval:0.01.02"
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

task Quast {
    input {
        Array[File] asm
        Array[String] asm_name

        File ref_fasta
        File gff

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(asm, "GiB") + size(ref_fasta, "GiB") + size(gff, "GiB"))

    command <<<
        set -euxo pipefail

        mkdir quast_results
        quast -o ./quast_results -r ~{ref_fasta} -g ~{gff} -l "~{sep=', ' asm_name}" -e --circos ~{sep=' ' asm}

        tar -cvzf quast_results.tar.gz quast_results
    >>>

    output {
        File report_html           = "quast_results/report.html"
        File report_tsv            = "quast_results/report.tsv"
        File report_txt            = "quast_results/report.txt"
        File report_pdf            = "quast_results/report.pdf"
        File transposed_report_tsv = "quast_results/transposed_report.tsv"
        File transposed_report_txt = "quast_results/transposed_report.txt"
        File quast_results         = "quast_results.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm-eval:0.01.02"
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


task PlotAsmVsRef {
    input {
        File asm
        String asm_name

        File ref_fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(asm, "GiB") + size(ref_fasta, "GiB"))

    command <<<
        set -euxo pipefail

        nucmer ~{ref_fasta} ~{asm} --delta=~{asm_name}.delta
        mummerplot --color --medium --filter --layout -R ~{ref_fasta} -Q ~{asm} --prefix ~{asm_name} --fat --png ~{asm_name}.delta

        df -h .
        du -hcs .
        find . type -f -exec ls -lah {} \;
    >>>

    output {
        File delta = "~{asm_name}.delta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm-eval:0.01.02"
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
