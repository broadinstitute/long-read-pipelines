version 1.0

import "../../../tasks/Phasing/WhatsHap.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"

workflow Statistics {
    meta{
        description : "..."
    }
    parameter_meta {
        one_chr_bams_from_all_samples:  "GCS path to subset BAM files"
        one_chr_bais_from_all_samples:  "GCS path to subset BAI file indices"
    }

    input {
        File phased_scaffold
        File phased_scaffold_index
        String chromosome
        Array[File] SampleIDs
    }
    
    scatter (sampleid in SampleIDs) {
        call SplitJointCallbySample.SplitVCFbySample as SP { input:
            joint_vcf = phased_scaffold,
            joint_vcf_tbi = phased_scaffold_index,
            region = chromosome,
            samplename = sampleid
        }

        call WhatsHap.Stats as WS { input:
            phased_vcf = SP.single_sample_vcf,
            phased_tbi = SP.single_sample_vcf_tbi
        }
    }

    output{
        Array[File] whatshap_statistics = WS.stats_tsv
    }
}
