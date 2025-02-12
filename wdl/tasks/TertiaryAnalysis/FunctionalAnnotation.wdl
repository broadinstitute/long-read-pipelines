version 1.0

import "../../structs/Structs.wdl"

task FunctionallyAnnotateVariants {
    input {
        File vcf
        File snpeff_db

        String snpeff_db_identifier

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 5*ceil(size([vcf, snpeff_db], "GB"))
    String prefix = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        set -x

        gunzip -c ~{snpeff_db} | tar xvf -

        /snpEff/scripts/snpEff ann -v \
            -c $PWD/snpeff_db/snpEff.config \
            -dataDir $PWD/snpeff_db/data \
            ~{snpeff_db_identifier} \
            ~{vcf} | bgzip > ~{prefix}.annotated.vcf.gz

        mv snpEff_summary.html ~{prefix}.snpEff_summary.html
        mv snpEff_genes.txt ~{prefix}.snpEff_genes.txt

        # Index the output VCF file:
        tabix -p vcf ~{prefix}.annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_index = "~{prefix}.annotated.vcf.gz.tbi"
        File snpEff_summary = "~{prefix}.snpEff_summary.html"
        File snpEff_genes = "~{prefix}.snpEff_genes.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-functional-annotation:0.0.1"
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
