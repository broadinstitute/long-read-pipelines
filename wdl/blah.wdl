version 1.0

import "tasks/Utils.wdl"
import "tasks/Finalize.wdl" as FF
import "tasks/utils/CloudUtils.wdl"
import "tasks/VariantUtils.wdl"

workflow blah {
    input {
        File VCF
        File idx
        String new_sample_name
    }

    call VariantUtils.ReSampleSingleSampleVCF as ReSample {input: vcf = VCF, idx = idx, new_sample_name = new_sample_name, out_prefix = "does_not_matter"}
    call CloudUtils.GetBlobFolder {input: blob = VCF}
    String outfolder = GetBlobFolder.bucket_and_folder
    String final_out_prefix = basename(basename(VCF, '.vcf.gz'), ".vcf")
    call FF.FinalizeToFile as ff_vcf {input: file = ReSample.reheadered_vcf, outdir = outfolder, name = final_out_prefix + ".RH.vcf.gz"}
    call FF.FinalizeToFile as ff_tbi {input: file = ReSample.reheadered_tbi, outdir = outfolder, name = final_out_prefix + ".RH.vcf.gz.tbi"}

    output {
        File rh_vcf = ff_vcf.gcs_path
        File rh_tbi = ff_tbi.gcs_path
    }
}