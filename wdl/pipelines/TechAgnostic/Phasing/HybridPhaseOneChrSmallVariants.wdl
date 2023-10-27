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
        use_margin: "Whatshap, Margin"
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
        Boolean use_margin = true
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
        
        call Utils.InferSampleName as InferSampleNamewhatshap { input: 
            bam = bam, 
            bai = bai}

        String sample_id = InferSampleNamewhatshap.sample_name

        call SplitJointCallbySample.SplitVCFbySample as Splitwhatshap { input:
            joint_vcf = joint_call.joint_gvcf,
            region = chromosome,
            samplename = sample_id
        }

        ########### call whatshap phasing the second time
        if (use_margin) {
            call MarginPhase.MarginPhase as M { input:
                bam = bam,
                bai = bai,
                data_type = "PacBio",
                unphased_vcf = Splitwhatshap.single_sample_vcf,
                unphased_tbi = Splitwhatshap.single_sample_vcf_tbi,
                ref_fasta = reference,
                ref_fasta_fai = reference_index
            }
        }
        if (!use_margin) {
            call WhatsHap.Phase as W { input:
                bam = bam,
                bai = bai,
                ref_fasta = reference,
                ref_fasta_fai = reference_index,
                unphased_vcf = Splitwhatshap.single_sample_vcf,
                unphased_tbi = Splitwhatshap.single_sample_vcf_tbi,
                chromosome = chromosome
            }
        }
        ########### call margin phasing the second time
       
        ######## Done calling physical phasing #########
    }


    # !CHOICE! 
    #############
    Array[File?] phased_snp_vcf = if (use_margin) then W.phased_vcf else M.phased_vcf
    Array[File?] phased_snp_tbi = if (use_margin) then W.phased_tbi else M.phased_tbi


    ##########Merge vcfs given different output##########   
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSampleswhatshap { input:
        vcf_input = phased_snp_vcf,
        tbi_input = phased_snp_tbi,
        pref = prefix
    }
    ################ Done Merging #################
    call StatPhase.Shapeit4 as Shapeit4 { input:
        vcf_input = MergeAcrossSampleswhatshap.merged_vcf,
        vcf_index = MergeAcrossSampleswhatshap.merged_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome,
        num_threads = num_t
    }
    ############# Done Statistical Phasing ################

    ############## Finalize##############
    call FF.FinalizeToFile as FinalizeGVCFwhatshap { input: outdir = gcs_output_directory, file = Shapeit4.scaffold_vcf }
    
    output{
        File phased_whatshap_scaffold = Shapeit4.scaffold_vcf
    }
}
