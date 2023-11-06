version 1.0

import "../../../tasks/Phasing/WhatsHap.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"
import "../../../tasks/Utility/VariantUtils.wdl" as VU

workflow Statistics {
    meta{
        description : "..."
    }
    parameter_meta {
        
    }

    input {
        File phased_vcf
        File phased_vcf_tbi
        String chromosome
    }
    
    call VU.SubsetVCF as SP { input:
        vcf_gz = phased_vcf,
        vcf_tbi = phased_vcf_tbi,
        locus = chromosome
    }

    call WhatsHap.Stats as WS { input:
        phased_vcf = SP.subset_vcf,
        phased_tbi = SP.subset_tbi
    }

    output{
        File whatshap_statistics = WS.stats_tsv
    }
}
