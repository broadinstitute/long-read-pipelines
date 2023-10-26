version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/WhatsHap.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"
import "../VariantCalling/LRJointCallGVCFs.wdl" as Jointcall
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Phasing/MarginPhase.wdl"
import "../../../tasks/Phasing/Hiphase.wdl"


workflow HybridPhase {
    meta{
        description : "..."
    }
    parameter_meta {
        one_chr_bams_from_all_samples:  "GCS path to subset BAM files"
        one_chr_bais_from_all_samples:  "GCS path to subset BAI file indices"
        one_chr_joint_vcf:  "path to subset joint vcf per chromosome"
        one_chr_joint_vcf_tbi:  "path to subset joint vcf index per chromosome"
        reference: "path to reference genome fasta file"
        reference_index: "path to reference genome index fai file"
        genetic_mapping_tsv_for_shapeit4: "path to the tsv file for the genetic mapping file address per chromosome"
        chromosome: "string for chromosome to be processed"
        prefix: "output file prefix, usually using chromosome"
        physical_phasing_type: "Whatshap, Margin, Hiphase"
        num_t: "integer for threads"
    }

    input {
        Array[File] one_chr_bams_from_all_samples
        Array[File] one_chr_bais_from_all_samples
        Array[File] one_chr_g_dv_vcf
        Array[File] one_chr_g_dv_vcf_tbi
        File reference
        File reference_index
        File ref_map
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String prefix
        String gcs_output_directory
        String physical_phasing_type
        Int num_t
    }
    
    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    # jointcall
    call Jointcall.LRJointCallGVCFs as joint_call { input:
        gvcfs = one_chr_g_dv_vcf,
        tbis = one_chr_g_dv_vcf_tbi,
        ref_map_file = ref_map, 
        prefix = chromosome,
        gcs_out_root_dir = gcs_output_directory
    }

    scatter (bam_bai in zip(one_chr_bams_from_all_samples, one_chr_bais_from_all_samples)) {
        File bam = bam_bai.left
        File bai = bam_bai.right
        
        call Utils.InferSampleName { input: 
            bam = bam, 
            bai = bai}
        String sample_id = InferSampleName.sample_name

        call SplitJointCallbySample.SplitVCFbySample as Split { input:
            joint_vcf = joint_call.joint_gvcf,
            region = chromosome,
            samplename = sample_id
        }

        if (physical_phasing_type == "Whatshap") {
            call WhatsHap.Phase as W { input:
                bam = bam,
                bai = bai,
                ref_fasta = reference,
                ref_fasta_fai = reference_index,
                unphased_vcf = Split.single_sample_vcf,
                unphased_tbi = Split.single_sample_vcf_tbi,
                chromosome = chromosome
                }
        }######## Done calling physical phasing #########

        if (physical_phasing_type == "Margin") {
            call MarginPhase.MarginPhase as M { input:
                bam = bam,
                bai = bai,
                data_type = "PacBio",
                unphased_vcf = Split.single_sample_vcf,
                unphased_tbi = Split.single_sample_vcf_tbi,
                ref_fasta = reference,
                ref_fasta_fai = reference_index
            }
        }######## Done calling physical phasing #########


    }
    ##########Merge vcfs given different output##########   
    if (physical_phasing_type == "Whatshap") {
        call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
                vcf_input = W.phased_vcf,
                tbi_input = W.phased_tbi,
                pref = prefix
        }
    }
    if (physical_phasing_type == "Margin") {
        call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
                vcf_input = M.phased_vcf,
                tbi_input = M.phased_tbi,
                pref = prefix
        }        
    }
    ################ Done Merging #################

    ############# Statistical Phasing ################
    call StatPhase.Shapeit4 as Shapeit4 { input:
        vcf_input = MergeAcrossSamples.merged_vcf,
        vcf_index = MergeAcrossSamples.merged_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome,
        num_threads = num_t
    }

    ############## Finalize##############
    call FF.FinalizeToFile as FinalizeGVCF { input: outdir = gcs_output_directory, file = Shapeit4.scaffold_vcf }
    
    output{
        File phased_scaffold = Shapeit4.scaffold_vcf
    }
}
