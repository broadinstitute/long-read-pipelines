version 1.0

import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Utility/VariantUtils.wdl" as VU

workflow Shapeit4PhaseWholeGenome {
    meta{
        description: "a workflow that extract vcfs and bams in a genomic interval"
    }
    input{
        Array[String] chr_list
        File joint_vcf
        File joint_vcf_tbi
        File genetic_mapping_tsv
        Int num_t
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv)

    # double scatter: first by chr, then by sample
    scatter (genome_region in chr_list) {
        call VU.SubsetVCF as SubsetVCF { input:
            vcf_gz = joint_vcf,
            vcf_tbi = joint_vcf_tbi,
            locus = genome_region
        }

        call StatPhase.Shapeit4 { input:
            vcf_input = SubsetVCF.subset_vcf,
            vcf_index = SubsetVCF.subset_tbi,
            mappingfile = genetic_mapping_dict[genome_region],
            region = genome_region,
            num_threads = num_t
        }

    }

   

    output{
        Array[File] phased_scaffold = Shapeit4.scaffold_vcf
    }


}


