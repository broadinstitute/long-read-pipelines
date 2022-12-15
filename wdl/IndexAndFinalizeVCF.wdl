version 1.0

import "tasks/Finalize.wdl" as FF
import "tasks/VariantUtils.wdl" as VU

workflow IndexAndFinalizeVCF {
    input {
        File vcf
        String sample_name
        String gcs_out_root_dir
    }

    call VU.IndexVCF { input: vcf = vcf }

    String dir = sub(gcs_out_root_dir, "/$", "") + "/PAV/~{sample_name}"
    call FF.FinalizeToFile as FinalizeVCF { input: outdir = dir, file = vcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = dir, file = IndexVCF.tbi }
    output {
        File finalized_vcf = FinalizeVCF.gcs_path
        File finalized_tbi = FinalizeTBI.gcs_path
    }
}
