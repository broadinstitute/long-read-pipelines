version 1.0

import "../../../tasks/VariantCalling/DeepVariant.wdl" as DeepVariant

workflow DeepTrio {
    input {
        File parent1_bam
        File parent1_bai
        String parent1_sample_id

        File parent2_bam
        File parent2_bai
        String parent2_sample_id

        File proband_bam
        File proband_bai
        String proband_sample_id

        File ref_fasta
        File ref_fasta_fai
        File ref_fasta_dict

        String model_type = "WGS"
    }

    call DeepVariant.DeepTrio as t_001_DeepTrio {
        input:
            parent1_bam = parent1_bam,
            parent1_bai = parent1_bai,
            parent1_sample_id = parent1_sample_id,

            parent2_bam = parent2_bam,
            parent2_bai = parent2_bai,
            parent2_sample_id = parent2_sample_id,

            proband_bam = proband_bam,
            proband_bai = proband_bai,
            proband_sample_id = proband_sample_id,

            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,

            model_type = model_type,
    }

    output {
        File parent1_vcf = t_001_DeepTrio.parent1_vcf
        File parent1_vcf_index = t_001_DeepTrio.parent1_vcf_index
        File parent1_gvcf = t_001_DeepTrio.parent1_gvcf
        File parent1_gvcf_index = t_001_DeepTrio.parent1_gvcf_index

        File parent2_vcf = t_001_DeepTrio.parent2_vcf
        File parent2_vcf_index = t_001_DeepTrio.parent2_vcf_index
        File parent2_gvcf = t_001_DeepTrio.parent2_gvcf
        File parent2_gvcf_index = t_001_DeepTrio.parent2_gvcf_index

        File proband_vcf = t_001_DeepTrio.proband_vcf
        File proband_vcf_index = t_001_DeepTrio.proband_vcf_index
        File proband_gvcf = t_001_DeepTrio.proband_gvcf
        File proband_gvcf_index = t_001_DeepTrio.proband_gvcf_index

        File logs = t_001_DeepTrio.logs
    }
}