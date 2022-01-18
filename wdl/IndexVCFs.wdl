version 1.0

import "tasks/VariantUtils.wdl"

import "tasks/Finalize.wdl" as FF

workflow IndexVCFs {
    input {
        File sniffles_vcf
        File pbsv_vcf

        String out_root
        String wgs_workflow_name
        String sample_name
    }
    String dir = out_root + '/~{wgs_workflow_name}/~{sample_name}/variants/sv/'

    call VariantUtils.IndexVCF as IndexSniffles { input: vcf = sniffles_vcf}
    call VariantUtils.IndexVCF as IndexPBSV     { input: vcf = pbsv_vcf}

    call FF.FinalizeToFile as FinalizeSniffles { input: outdir = dir, file = IndexSniffles.tbi }
    call FF.FinalizeToFile as FinalizePBSV     { input: outdir = dir, file = IndexPBSV.tbi }

    output {
        File sniffles_tbi = FinalizeSniffles.gcs_path
        File pbsv_tbi = FinalizePBSV.gcs_path
    }
}