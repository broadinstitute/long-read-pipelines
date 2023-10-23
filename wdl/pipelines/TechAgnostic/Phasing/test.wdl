version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/MarginPhase.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"

workflow test{
    input {
        Array[File] VCFs
        Array[File] TBIs
        String pref

    }

    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
        vcf_input = VCFs,
        tbi_input = TBIs,       
        pref = pref
    }
    output{
        File merged = MergeAcrossSamples.merged_vcf
    }
}