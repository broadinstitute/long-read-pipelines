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
        File wholegenome_bams
        File wholegenome_bais
        File wholegenome_joint_vcf
        File wholegenome_joint_vcf_tbi
        File one_chr_joint_sv
        File one_chr_joint_sv_tbi
        File reference
        File reference_index
        String chromosome
        String prefix
        Int num_t
    }


    call VU.SubsetVCF as SubsetSNPsJoint { input:
            vcf_gz = wholegenome_joint_vcf,
            vcf_tbi = wholegenome_joint_vcf_tbi,
            locus = chromosome
        }

    # call VU.SubsetVCF as SubsetSVsJoint { input:
    #     vcf_gz = wholegenome_joint_sv,
    #     vcf_tbi = wholegenome_joint_sv_tbi,
    #     locus = chromosome
    # }

 
    ####### Subset Bam####
    call U.SubsetBam as SubsetBam { input:
        bam = wholegenome_bams,
        bai = wholegenome_bais,
        locus = chromosome
    }
    ####### DONE ######


    call U.InferSampleName { input: 
        bam = wholegenome_bams, 
        bai = wholegenome_bais
    }
    String sample_id = InferSampleName.sample_name

    call SplitJointCallbySample.SplitVCFbySample as SP { input:
        joint_vcf = SubsetSNPsJoint.subset_vcf,
        region = chromosome,
        samplename = sample_id
    }

    call SplitJointCallbySample.SplitVCFbySample as SV_split { input:
        joint_vcf = one_chr_joint_sv,
        region = chromosome,
        samplename = sample_id
    }

    call Hiphase.HiphaseSVs as HP_SV { input:
        bam = SubsetBam.subset_bam,
        bai = SubsetBam.subset_bai,
        unphased_snp_vcf = SP.single_sample_vcf,
        unphased_snp_tbi = SP.single_sample_vcf_tbi,
        unphased_sv_vcf = SV_split.single_sample_vcf,
        unphased_sv_tbi = SV_split.single_sample_vcf_tbi,
        ref_fasta = reference,
        ref_fasta_fai = reference_index,
        samplename = sample_id

    }
    

    output{
        File hiphase_snp = HP_SV.phased_snp_vcf
        File hiphase_snp_tbi = HP_SV.phased_snp_vcf_tbi
        File hiphase_sv = HP_SV.phased_sv_vcf
        File hiphase_sv_tbi = HP_SV.phased_sv_vcf_tbi
        

    }
}
