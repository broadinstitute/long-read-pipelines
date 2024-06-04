version 1.0

import "HybridPhaseOneChrSmallVariantHyphase.wdl" as PhaseOneChr

workflow HybridPhaseWholeGenome {
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        String output_dir
        Array[String] chr_list
        Array[File] whole_genome_bams
        Array[File] whole_genome_bais
        File joint_vcf
        File joint_vcf_tbi
        File reference
        File reference_index
        File genetic_mapping_tsv
        Int num_t
    }

    # Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv)

    # double scatter: first by chr, then by sample
    scatter (genome_region in chr_list) {
        call PhaseOneChr.HybridPhase as HybridPhase { input:
            gcs_out_root_dir = output_dir,
            wholegenome_bams_from_all_samples = whole_genome_bams,
            wholegenome_bais_from_all_samples = whole_genome_bais,
            wholegenome_joint_vcf = joint_vcf,
            wholegenome_joint_vcf_tbi = joint_vcf_tbi,
            reference = reference,
            reference_index = reference_index,
            genetic_mapping_tsv_for_shapeit4 = genetic_mapping_tsv,
            chromosome = genome_region,
            prefix = genome_region,
            num_t = num_t

        }

    }   

    output{
        Array[File] phased_snp_scaffold = HybridPhase.snp_shapeit4scaffold
        Array[File] phased_sv_svaffold = HybridPhase.sv_shapeit4scaffold
    }


}


