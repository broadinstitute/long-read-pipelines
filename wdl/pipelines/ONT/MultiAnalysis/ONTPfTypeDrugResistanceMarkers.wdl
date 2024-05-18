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
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."

        do_functional_annotation: "Whether to perform functional annotation"
    }

    input {
        File vcf
        File snpeff_db
        File drug_resistance_list

        String dir_prefix
        String? gcs_out_root_dir

        Boolean do_functional_annotation = true
    }

    if (do_functional_annotation) {
        call FUNK.FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db }
    }

    call CallDrugResistanceMutations {
        input:
            vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf]),
            drug_resistance_list = drug_resistance_list
    }

    # Finalize data
    if (defined(gcs_out_root_dir)) {

        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])

        String outdir = sub(concrete_gcs_out_root_dir, "/$", "") + "/ONTPfTypeDrugResistanceMarkers/~{dir_prefix}"
        String dir = outdir + "/reports"

        call FF.FinalizeToFile as FinalizeDRReport { input: outdir = dir, file = CallDrugResistanceMutations.report }

        if (do_functional_annotation) {
            call FF.FinalizeToFile as FinalizeAnnotatedVCF { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf]) }
            call FF.FinalizeToFile as FinalizeAnnotatedVCFIndex { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf_index]) }
            call FF.FinalizeToFile as FinalizeSnpEffSummary { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_summary]) }
            call FF.FinalizeToFile as FinalizeSnpEffGenes { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_genes]) }
        }

        File final_annotated_vcf = select_first([FinalizeAnnotatedVCF.gcs_path, FunctionallyAnnotateVariants.annotated_vcf])
        File final_annotated_vcf_index = select_first([FinalizeAnnotatedVCFIndex.gcs_path, FunctionallyAnnotateVariants.annotated_vcf_index])
        File final_snpeff_summary = select_first([FinalizeSnpEffSummary.gcs_path, FunctionallyAnnotateVariants.snpEff_summary])
        File final_snpeff_genes = select_first([FinalizeSnpEffGenes.gcs_path, FunctionallyAnnotateVariants.snpEff_genes])
    }

    output {
        File drug_res_report = select_first([FinalizeDRReport.gcs_path, CallDrugResistanceMutations.report])

        File? annotated_vcf = final_annotated_vcf
        File? annotated_vcf_index = final_annotated_vcf_index
        File? snpEff_summary = final_snpeff_summary
        File? snpEff_genes = final_snpeff_genes
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
