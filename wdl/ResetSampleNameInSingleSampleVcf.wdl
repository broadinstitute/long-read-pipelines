version 1.0

import "tasks/Utils.wdl"
import "tasks/Finalize.wdl" as FF
import "tasks/utils/CloudUtils.wdl"
import "tasks/VariantUtils.wdl"

import "blah.wdl"

workflow ResetSampleNameInSingleSampleVcf {
    meta {
        description: "For resetting the sample name in a single-sample VCF file to the requested value"
    }

    input {
        File dvp_g_tbi
        File dvp_g_vcf
        File dvp_phased_tbi
        File dvp_phased_vcf
        File dvp_tbi
        File dvp_vcf
        File pbsv_tbi
        File pbsv_vcf
        File sniffles_tbi
        File sniffles_vcf

        String new_sample_name
    }

    parameter_meta {
        VCF: "VCF to operate on."
        idx: "Accompanying index for VCF."
        new_sample_name: "New sample name desired (assumed to be a valid value)."
    }

    call blah.blah as DVP_g      {input: VCF = dvp_g_vcf, idx = dvp_g_tbi, new_sample_name = new_sample_name}
    call blah.blah as DVP_phased {input: VCF = dvp_phased_vcf, idx = dvp_phased_tbi, new_sample_name = new_sample_name}
    call blah.blah as DVP        {input: VCF = dvp_vcf, idx = dvp_tbi, new_sample_name = new_sample_name}
    call blah.blah as PBSV       {input: VCF = pbsv_vcf, idx = pbsv_tbi, new_sample_name = new_sample_name}
    call blah.blah as Sniffles   {input: VCF = sniffles_vcf, idx = sniffles_tbi, new_sample_name = new_sample_name}

    output {
        File rh_dvp_g_tbi = DVP_g.rh_tbi
        File rh_dvp_g_vcf = DVP_g.rh_vcf
        File rh_dvp_phased_tbi = DVP_phased.rh_tbi
        File rh_dvp_phased_vcf = DVP_phased.rh_vcf
        File rh_dvp_tbi = DVP.rh_tbi
        File rh_dvp_vcf = DVP.rh_vcf
        File rh_pbsv_tbi = PBSV.rh_tbi
        File rh_pbsv_vcf = PBSV.rh_vcf
        File rh_sniffles_tbi = Sniffles.rh_tbi
        File rh_sniffles_vcf = Sniffles.rh_vcf

        # purely Terra trick
        File orig_dvp_g_tbi = dvp_g_tbi
        File orig_dvp_g_vcf = dvp_g_vcf
        File orig_dvp_phased_tbi = dvp_phased_tbi
        File orig_dvp_phased_vcf = dvp_phased_vcf
        File orig_dvp_tbi = dvp_tbi
        File orig_dvp_vcf = dvp_vcf
        File orig_pbsv_tbi = pbsv_tbi
        File orig_pbsv_vcf = pbsv_vcf
        File orig_sniffles_tbi = sniffles_tbi
        File orig_sniffles_vcf = sniffles_vcf
    }
}
