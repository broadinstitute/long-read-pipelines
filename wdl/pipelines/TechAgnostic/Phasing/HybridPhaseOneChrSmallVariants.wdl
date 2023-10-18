version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/WhatsHap.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"


workflow HybridPhase {
    meta{
        description : "..."
    }
    parameter_meta {
        one_chr_bams_from_all_samples:  "GCS path to subset BAM files"
        one_chr_bais_from_all_samples:  "GCS path to subset BAI file indices"
    }

    input {
        Array[File] one_chr_bams_from_all_samples
        Array[File] one_chr_bais_from_all_samples
        File one_chr_joint_vcf
        File one_chr_joint_vcf_tbi
        File reference
        File reference_index
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String prefix
        Int num_t
    }
    
    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    scatter (bam_bai in zip(one_chr_bams_from_all_samples, one_chr_bais_from_all_samples)) {
        File bam = bam_bai.left
        File bai = bam_bai.right
        
        call Utils.InferSampleName { input: 
            bam = bam, 
            bai = bai}
        String sample_id = InferSampleName.sample_name

        call SplitJointCallbySample.SplitVCFbySample as Split { input:
            joint_vcf = one_chr_joint_vcf,
            region = chromosome,
            samplename = sample_id
        }

        call WhatsHap.Phase as whatshap_phasing { input:
            bam = bam,
            bai = bai,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            unphased_vcf = Split.single_sample_vcf,
            unphased_tbi = Split.single_sample_vcf_tbi,
            chromosome = chromosome
        }
    }
        
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
        vcf_input = whatshap_phasing.phased_vcf,
        tbi_input = whatshap_phasing.phased_tbi,
        pref = prefix
    }

    call StatPhase.Shapeit4 { input:
        vcf_input = MergeAcrossSamples.merged_vcf,
        vcf_index = MergeAcrossSamples.merged_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome,
        num_threads = num_t
    }

    output{
        File phased_scaffold = Shapeit4.scaffold_vcf
    }
}
