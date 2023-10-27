version 1.0

import "../../../tasks/Utility/Utils.wdl" as U
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/Hiphase.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF


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
        num_t: "integer for threads"
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
        String gcs_out_root_dir
        Int num_t
    }
    
    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    scatter (bam_bai in zip(one_chr_bams_from_all_samples, one_chr_bais_from_all_samples)) {
        File bam = bam_bai.left
        File bai = bam_bai.right
        
        call U.InferSampleName { input: 
            bam = bam, 
            bai = bai}
        String sample_id = InferSampleName.sample_name

        call SplitJointCallbySample.SplitVCFbySample as SP { input:
            joint_vcf = one_chr_joint_vcf,
            region = chromosome,
            samplename = sample_id
        }

        call Hiphase.HiphaseSNPs as HP { input:
            bam = bam,
            bai = bai,
            unphased_snp_vcf = SP.single_sample_vcf,
            unphased_snp_tbi = SP.single_sample_vcf_tbi,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            samplename = sample_id

        }
    }
        
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
        vcf_input = HP.phased_vcf,
        tbi_input = HP.phased_vcf_tbi,
        pref = prefix
    }

    # call StatPhase.Shapeit4 { input:
    #     vcf_input = MergeAcrossSamples.merged_vcf,
    #     vcf_index = MergeAcrossSamples.merged_tbi,
    #     mappingfile = genetic_mapping_dict[chromosome],
    #     region = chromosome,
    #     num_threads = num_t
    # }
    call FF.FinalizeToDir as Finalizevcfs {
        input: outdir = gcs_out_root_dir, files = HP.phased_vcf
    }
    call FF.FinalizeToDir as Finalizetbis {
        input: outdir = gcs_out_root_dir, files = HP.phased_vcf_tbi
    }

    output{
        Array[File] hiphased_vcf = HP.phased_vcf
        Array[File] hiphased_vcf_tbi = HP.phased_vcf_tbi
        File merged_vcf = MergeAcrossSamples.merged_vcf
    }
}
