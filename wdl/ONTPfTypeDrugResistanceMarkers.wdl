version 1.0

import "tasks/Structs.wdl"
import "tasks/Finalize.wdl" as FF

workflow ONTPfTypeDrugResistanceMarkers {
    input {
        File vcf
        File snpeff_db
        File drug_resistance_list

        String dir_prefix
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTPfTypeDrugResistanceMarkers/~{dir_prefix}"

    call FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db }

    call CallDrugResistanceMutations {
        input:
            vcf = FunctionallyAnnotateVariants.annotated_vcf,
            drug_resistance_list = drug_resistance_list
    }

    # Finalize data
    String dir = outdir + "/reports"

    call FF.FinalizeToFile as FinalizeAnnotatedVCF { input: outdir = dir, file = FunctionallyAnnotateVariants.annotated_vcf }
    call FF.FinalizeToFile as FinalizeDRReport { input: outdir = dir, file = CallDrugResistanceMutations.report }

    output {
        File annotated_vcf = FinalizeAnnotatedVCF.gcs_path
        File drug_res_report = FinalizeDRReport.gcs_path
    }
}

task FunctionallyAnnotateVariants {
    input {
        File vcf
        File snpeff_db

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([vcf, snpeff_db], "GB"))
    String prefix = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        set -x

        gunzip -c ~{snpeff_db} | tar xvf -

        /usr/local/bin/snpEff ann -v \
            -c $PWD/snpeff_db/snpEff.config \
            -dataDir $PWD/snpeff_db/data \
            PlasmoDB-61_Pfalciparum3D7_Genome \
            ~{vcf} \
            | gzip > ~{prefix}.annotated.vcf.gz

        mv snpEff_summary.html ~{prefix}.snpEff_summary.html
        mv snpEff_genes.txt ~{prefix}.snpEff_genes.txt
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File snpEff_summary = "~{prefix}.snpEff_summary.html"
        File snpEff_genes = "~{prefix}.snpEff_genes.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/snpeff:5.1d--hdfd78af_0"
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

task CallDrugResistanceMutations {
    input {
        File vcf
        File drug_resistance_list

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([vcf, drug_resistance_list], "GB"))
    String prefix = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        set -x

        while read LINE; do
            GENE_NAME=$(echo $LINE | awk '{print $1}')
            GENE_ID=$(echo $LINE | awk '{print $2}')
            MUTATION=$(echo $LINE | awk '{print $3}')

            grep $GENE_ID ~{vcf} | \
                grep $MUTATION | \
                wc -l | \
                awk -v gene_name=$GENE_NAME \
                    -v gene_id=$GENE_ID \
                    -v mutation=$MUTATION \
                    '{ gene_name "\t" gene_id "\t" mutation "\t" ($1 > 0) }' \
            >> ~{prefix}.drug_resistance_report.txt
        done <~{drug_resistance_list}
    >>>

    output {
        File report = "~{prefix}.drug_resistance_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/snpeff:5.1d--hdfd78af_0"
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
