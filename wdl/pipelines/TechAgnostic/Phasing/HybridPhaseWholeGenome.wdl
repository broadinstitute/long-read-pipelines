version 1.0

import "HybridPhaseOneChrSmallVariants.wdl" as PhaseOneChr
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Utility/Utils.wdl" as U

workflow HybridPhaseWholeGenome {
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        Array[File] whole_genome_bams
        Array[File] whole_genome_bais
        Array[String] sampleIds
        File joint_vcf
        File joint_vcf_tbi
        File reference
        File reference_index
        File genetic_mapping_tsv
        String prefix
        Int num_t
    }
    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv)

    Array[String] chr_list = ["chr1", "chr6",]

    # double scatter: first by chr, then by sample
    scatter (genome_region in chr_list) {
        call VU.SubsetVCF as SubsetVCF { input:
            vcf_gz = joint_vcf,
            vcf_tbi = joint_vcf_tbi,
            locus = genome_region
        }

        scatter (bam_bai in zip(whole_genome_bams, whole_genome_bais)) {
            File bam = bam_bai.left
            File bai = bam_bai.right

            call U.SubsetBam as SubsetBam { input:
                bam = bam,
                bai = bai,
                locus = genome_region
            }
        }   

        call PhaseOneChr.HybridPhase as HybridPhase { input:
            one_chr_bams_from_all_samples = SubsetBam.subset_bam,
            one_chr_bais_from_all_samples = SubsetBam.subset_bai,
            one_chr_joint_vcf = SubsetVCF.subset_vcf,
            one_chr_joint_vcf_tbi = SubsetVCF.subset_tbi,
            reference = reference,
            reference_index = reference_index,
            genetic_mapping_tsv_for_shapeit4 = genetic_mapping_tsv,
            chromosome = genome_region,
            prefix = genome_region,
            num_t = num_t
        }

    }   

    output{
        Array[File] phased_scaffold = HybridPhase.phased_scaffold
    }


}


