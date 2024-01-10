version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/TertiaryAnalysis/FunctionalAnnotation.wdl" as FUNK
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ONTPfTypeDrugResistanceMarkers {

    meta {
        description: "Workflow to generate a report of drug resistance markers"
    }
    parameter_meta {
        vcf: "VCF file to process"
        snpeff_db: "SnpEff database for functional annotation"
        drug_resistance_list: "List of drug resistance markers for which to search"

        dir_prefix: "Prefix for output directory"
        gcs_out_root_dir: "GCS output root directory"

        do_functional_annotation: "Whether to perform functional annotation"
    }

    input {
        File vcf
        File snpeff_db
        File drug_resistance_list

        String dir_prefix
        String gcs_out_root_dir

        Boolean do_functional_annotation = true
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTPfTypeDrugResistanceMarkers/~{dir_prefix}"

    if (do_functional_annotation) {
        call FUNK.FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db }
    }

    call CallDrugResistanceMutations {
        input:
            vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf]),
            drug_resistance_list = drug_resistance_list
    }

    # Finalize data
    String dir = outdir + "/reports"

    call FF.FinalizeToFile as FinalizeDRReport { input: outdir = dir, file = CallDrugResistanceMutations.report }

    if (do_functional_annotation) {
        call FF.FinalizeToFile as FinalizeAnnotatedVCF { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf]) }
        call FF.FinalizeToFile as FinalizeAnnotatedVCFIndex { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf_index]) }
        call FF.FinalizeToFile as FinalizeSnpEffSummary { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_summary]) }
        call FF.FinalizeToFile as FinalizeSnpEffGenes { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_genes]) }
    }

    output {
        File drug_res_report = FinalizeDRReport.gcs_path

        File? annotated_vcf = FinalizeAnnotatedVCF.gcs_path
        File? annotated_vcf_index = FinalizeAnnotatedVCFIndex.gcs_path
        File? snpEff_summary = FinalizeSnpEffSummary.gcs_path
        File? snpEff_genes = FinalizeSnpEffGenes.gcs_path
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

            zcat ~{vcf} | grep $GENE_ID | grep $MUTATION | wc -l | \
                awk -v gene_name=$GENE_NAME -v gene_id=$GENE_ID -v mutation=$MUTATION \
                    '{ print gene_name, gene_id, mutation, ($1 > 0) ? "present" : "absent" }' | \
                tee -a ~{prefix}.drug_resistance_report.txt
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
