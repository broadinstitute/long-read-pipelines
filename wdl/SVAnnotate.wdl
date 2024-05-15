version 1.0

import "tasks/Structs.wdl"

workflow SVAnnotate {
    input {
        File vcf
        File gtf
        File ploidy
        File ref_fai
        String caller
        String prefix
    }

    parameter_meta {
        vcf: "VCF file to standardize"
        gtf: "Protein-coding gtf"
        ploidy: "ploidy table for sample, where rows are sample, chrs... and each line has ploidy number for each chr"
        ref_fai: "fai file for reference"
        caller: "caller used to generate VCF, for long-reads must be one of: pbsv, sniffles, pav, hapdiff"
        prefix: "Output file prefix"
    }

    call Standardize { input: vcf=vcf, ploidy=ploidy, prefix=prefix, ref_fai=ref_fai, caller=caller, gtf=gtf }

    output {
        File vcf_std = Standardize.vcf_std
        File vcf_std_tbi = Standardize.vcf_std_tbi
        File vcf_anno = Standardize.vcf_anno
        File vcf_anno_tbi = Standardize.vcf_anno_tbi
    }
}

task Standardize {
    input {
        File vcf
        File ploidy
        String prefix
        File ref_fai
        File gtf
        String caller

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 5 + 5*ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        tabix ~{vcf}

        svtk standardize --contigs ~{ref_fai} --min-size 25 --prefix ~{prefix} ~{vcf} - ~{caller} | bcftools sort /dev/stdin -o ~{prefix}.std.vcf.gz -O z
        tabix ~{prefix}.std.vcf.gz

        python /gatk/format_svtk_vcf_for_gatk.py --vcf ~{prefix}.std.vcf.gz --fix-end --out ~{prefix}.std.final.vcf.gz --ploidy-table ~{ploidy}
        tabix ~{prefix}.std.final.vcf.gz

        gatk SVAnnotate -V ~{prefix}.std.final.vcf.gz --protein-coding-gtf ~{gtf} -O ~{prefix}.anno.vcf
        bcftools view -Oz ~{prefix}.anno.vcf > ~{prefix}.anno.vcf.gz
        tabix ~{prefix}.anno.vcf.gz
    >>>

    output {
        File vcf_std = "~{prefix}.std.final.vcf.gz"
        File vcf_std_tbi = "~{prefix}.std.final.vcf.gz.tbi"
        File vcf_anno = "~{prefix}.anno.vcf.gz"
        File vcf_anno_tbi = "~{prefix}.anno.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-svannotate:latest"
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
