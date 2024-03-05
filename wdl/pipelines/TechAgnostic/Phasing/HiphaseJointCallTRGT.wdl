version 1.0

import "../../../tasks/Utility/Utils.wdl" as U
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/Hiphase.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"


workflow HiphaseJointCallTRGT {
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
        File trgt_vcf
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

    call index as ID { input:
        vcf = trgt_vcf,
        samplename = prefix
    }

    call VU.SubsetVCF as SubsetTrgtJoint { input:
        vcf_gz = trgt_vcf,
        vcf_tbi = ID.subset_tbi,
        locus = chromosome
    }

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

    call SplitJointCallbySample.SplitVCFbySample as Trgt_split { input:
        joint_vcf = trgt_vcf,
        region = chromosome,
        samplename = sample_id
    }

    call convert as cleanvcf{
        input:
            vcf = SV_split.single_sample_vcf,
            samplename = sample_id
    }

    call convert as cleantrgt{
        input:
            vcf = Trgt_split.single_sample_vcf,
            samplename = sample_id
    }

    call Hiphase.HiphaseSVTrgt as HP_SV_Trgt { input:
        bam = SubsetBam.subset_bam,
        bai = SubsetBam.subset_bai,
        unphased_snp_vcf = SP.single_sample_vcf,
        unphased_snp_tbi = SP.single_sample_vcf_tbi,
        unphased_sv_vcf = cleanvcf.subset_vcf,
        unphased_sv_tbi = cleanvcf.subset_tbi,
        unphased_trgt_vcf = cleantrgt.subset_vcf,
        unphased_trgt_tbi = cleantrgt.subset_tbi,
        ref_fasta = reference,
        ref_fasta_fai = reference_index,
        samplename = sample_id

    }
    

    output{
        File hiphase_snp = HP_SV_Trgt.phased_snp_vcf
        File hiphase_snp_tbi = HP_SV_Trgt.phased_snp_vcf_tbi
        File hiphase_sv = HP_SV_Trgt.phased_sv_vcf
        File hiphase_sv_tbi = HP_SV_Trgt.phased_sv_vcf_tbi
        File hiphase_trgt = HP_SV_Trgt.phased_trgt_vcf
        File hiphase_trgt_tbi = HP_SV_Trgt.phased_trgt_vcf_tbi
        File hiphase_haplotag = HP_SV_Trgt.haplotag_file
        

    }
}



task convert {

    meta {
        description: ""
    }

    parameter_meta {
        vcf: "VCF file to be cleaned"
    }

    input {
        File vcf
        String samplename
    }

    Int disk_size = 2*ceil(size([vcf], "GB")) + 1
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cp ~{docker_dir}/convert_lower_case.py ~{work_dir}/convert_lower_case.py
        cd ~{work_dir}

        python convert_lower_case.py -i ~{vcf} -o ~{samplename}_sv_cleaned.vcf
        bgzip ~{samplename}_sv_cleaned.vcf ~{samplename}_sv_cleaned.vcf.gz
        tabix -p vcf ~{samplename}_sv_cleaned.vcf.gz
    >>>

    output {
        File subset_vcf = "~{work_dir}/~{samplename}_sv_cleaned.vcf.gz"
        File subset_tbi = "~{work_dir}/~{samplename}_sv_cleaned.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 2
        memory:  "32 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"hangsuunc/cleanvcf:v1"
    }

}

task index{

    meta {
        description: ""
    }

    parameter_meta {
        vcf: "VCF file to be indexed"
    }

    input {
        File vcf
        String samplename
    }

    Int disk_size = 2*ceil(size([vcf], "GB")) + 1
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        
        tabix -p vcf ~{vcf}
    >>>

    output {
        File subset_tbi = "~{work_dir}/~{vcf}.tbi"
    }
    ###################
    runtime {
        cpu: 2
        memory:  "32 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"hangsuunc/cleanvcf:v1"
    }

}