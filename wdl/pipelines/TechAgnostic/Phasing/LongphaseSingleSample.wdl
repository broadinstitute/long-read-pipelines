version 1.0

import "../../../tasks/Utility/Utils.wdl" as U
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/Longphase.wdl"


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
        String samplename
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
    
    call Longphase.LongphaseSNPs as longphase { input:
            bam = SubsetBam.subset_bam,
            bai = SubsetBam.subset_bai,
            unphased_snp_vcf = SubsetSNPs.subset_vcf,
            unphased_snp_tbi = SubsetSNPs.subset_tbi,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            samplename = samplename
        }

    output{
        File phased_snp_vcf = longphase.phased_vcf
    }
}
