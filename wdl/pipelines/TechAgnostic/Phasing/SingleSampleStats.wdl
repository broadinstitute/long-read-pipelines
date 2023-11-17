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
        String chromosome
        String sampleID
    }

    call SplitJointCallbySample.SplitVCFbySample as SP { input:
            joint_vcf = phased_vcf,
            region = chromosome,
            samplename = sampleID
        }
    # call VU.SubsetVCF as SP { input:
    #     vcf_gz = Splitbysample.single_sample_vcf,
    #     vcf_tbi = Splitbysample.single_sample_vcf_tbi,
    #     locus = chromosome
    # }

    call WhatsHap.Stats as WS { input:
        phased_vcf = SP.single_sample_vcf,
        phased_tbi = SP.single_sample_vcf_tbi
    }

    output{
        File whatshap_statistics = WS.stats_tsv
    }
}
