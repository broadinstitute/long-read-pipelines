version 1.0

import "../../../tasks/Utility/Utils.wdl" as U
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/Hiphase.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"


workflow HybridPhase {
    meta{
        description : "..."
    }
    parameter_meta {
    }

    input {
        File all_chr_bams
        File all_chr_bais
        File deep_var_vcf
        File deep_var_vcf_tbi
        File pbsv_vcf
        File pbsv_vcf_tbi
        File reference
        File reference_index
        String chromosome
        String prefix
        Int num_t
    }

    call U.SubsetBam as SubsetBam { input:
            bam = all_chr_bams,
            bai = all_chr_bais,
            locus = chromosome
        }

    call VU.SubsetVCF as SubsetSNPs { input:
            vcf_gz = deep_var_vcf,
            vcf_tbi = deep_var_vcf_tbi,
            locus = chromosome
        }

    call VU.SubsetVCF as SubsetSVs { input:
        vcf_gz = pbsv_vcf,
        vcf_tbi = pbsv_vcf_tbi,
        locus = chromosome
    }
    
    call Hiphase.Hiphase as hiphase { input:
            bam = SubsetBam.subset_bam,
            bai = SubsetBam.subset_bai,
            unphased_snp_vcf = SubsetSNPs.subset_vcf,
            unphased_snp_tbi = SubsetSNPs.subset_tbi,
            unphased_sv_vcf = SubsetSVs.subset_vcf,
            unphased_sv_tbi = SubsetSVs.subset_tbi,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            prefix = prefix
        }

    output{
        File phased_snp_vcf = hiphase.phased_snp_vcf
        File phased_sv_vcf = hiphase.phased_sv_vcf
    }
}
