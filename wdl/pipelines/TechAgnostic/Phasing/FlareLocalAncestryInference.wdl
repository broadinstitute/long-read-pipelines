version 1.0

import "../../../tasks/Phasing/Flare.wdl"


workflow FlareLocalAncestryInference {
    meta{
        description : "..."
    }
    parameter_meta {
    }

    input {
        File ref_vcf
        File ref_vcf_index
        File ref_panel
        File test_vcf
        File test_vcf_index
        File plink_map
        String output_prefix
    }


    call Flare.Flare as F{
        input:
            ref_vcf = ref_vcf,
            ref_vcf_index = ref_vcf_index,
            ref_panel = ref_panel,
            test_vcf = test_vcf,
            test_vcf_index = test_vcf_index,
            plink_map = plink_map,
            output_prefix = output_prefix
    }

    output {
        Array[File] output_files = F.output_files
        

    }
}
