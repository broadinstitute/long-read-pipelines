version 1.0

import "tasks/Utils.wdl"
import "tasks/Finalize.wdl" as FF
import "tasks/utils/CloudUtils.wdl"
import "tasks/VariantUtils.wdl"

workflow ResetSampleNameInSingleSampleVcf {
    meta {
        description: "For resetting the sample name in a single-sample VCF file to the requested value"
    }

    input {
        File  VCF
        File? idx
        String new_sample_name
    }

    parameter_meta {
        VCF: "VCF to operate on."
        idx: "Accompanying index for VCF."
        new_sample_name: "New sample name desired (assumed to be a valid value)."
    }

    call VariantUtils.ReSampleSingleSampleVCF as ReSample {input: vcf = VCF, idx = idx, new_sample_name = new_sample_name, out_prefix = "does_not_matter"}

    call CloudUtils.GetBlobFolder {input: blob = VCF}
    String outfolder = GetBlobFolder.bucket_and_folder
    String final_out_prefix = basename(basename(VCF, '.vcf.gz'), ".vcf")
    call FF.FinalizeToFile as ff_vcf {input: file = ReSample.reheadered_vcf, outdir = outfolder, name = final_out_prefix + ".RH.vcf.gz"}
    call FF.FinalizeToFile as ff_tbi {input: file = ReSample.reheadered_tbi, outdir = outfolder, name = final_out_prefix + ".RH.vcf.gz.tbi"}

    output {
        File reheadered_vcf = ff_vcf.gcs_path
        File reheadered_tbi = ff_tbi.gcs_path

        # purely Terra trick
        File  original_vcf = VCF
        File? original_idx = idx
    }
}


