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
        Array[File] wholegenome_bams_from_all_samples
        Array[File] wholegenome_bais_from_all_samples

        File wholegenome_joint_vcf
        File wholegenome_joint_vcf_tbi
        File wholegenome_joint_sv
        File wholegenome_joint_sv_tbi
        File reference
        File reference_index
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String region
        String prefix
        String gcs_out_root_dir
        Int num_t
    }



    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)
    Int data_length = length(wholegenome_bams_from_all_samples)
    Array[Int] indexes= range(data_length)

    call VU.SubsetVCF as SubsetSNPsJoint { input:
            vcf_gz = wholegenome_joint_vcf,
            vcf_tbi = wholegenome_joint_vcf_tbi,
            locus = region
        }

    call VU.SubsetVCF as SubsetSVsJoint { input:
        vcf_gz = wholegenome_joint_sv,
        vcf_tbi = wholegenome_joint_sv_tbi,
        locus = region
    }

    scatter (idx in indexes)  {
        File all_chr_bam = wholegenome_bams_from_all_samples[idx]
        File all_chr_bai = wholegenome_bais_from_all_samples[idx]

        ####### Subset Bam####
        call U.SubsetBam as SubsetBam { input:
            bam = all_chr_bam,
            bai = all_chr_bai,
            locus = region
        }

        ####### DONE ######


        call U.InferSampleName { input: 
            bam = all_chr_bam, 
            bai = all_chr_bai
        }
        String sample_id = InferSampleName.sample_name

        call SplitJointCallbySample.SplitVCFbySample as SP { input:
            joint_vcf = SubsetSNPsJoint.subset_vcf,
            region = region,
            samplename = sample_id
        }

        call SplitJointCallbySample.SplitVCFbySample as SV_split { input:
            joint_vcf = SubsetSVsJoint.subset_vcf,
            region = region,
            samplename = sample_id
        }

        call convert as cleanvcf{
            input:
                vcf = SV_split.single_sample_vcf,
                samplename = sample_id
                
        }

        call Hiphase.HiphaseSVs as HP_SV { input:
            bam = SubsetBam.subset_bam,
            bai = SubsetBam.subset_bai,
            unphased_snp_vcf = SP.single_sample_vcf,
            unphased_snp_tbi = SP.single_sample_vcf_tbi,
            unphased_sv_vcf = cleanvcf.subset_vcf,
            unphased_sv_tbi = cleanvcf.subset_tbi,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            samplename = sample_id

        }
    }
    
    ### phase small variants as scaffold  
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
        vcf_input = HP_SV.phased_snp_vcf,
        tbi_input = HP_SV.phased_snp_vcf_tbi,
        pref = prefix
    }
    call SplitMultiallelicCalls as Norm_SNPs { input:
        bcftools_merged_vcf = MergeAcrossSamples.merged_vcf,
        bcftools_merged_vcf_tbi = MergeAcrossSamples.merged_tbi,
        prefix = prefix
    }
    ##### add filtering small variants step #######
    call FilterSmallVariants as Filter_SNPs{ input:
        bcftools_vcf = Norm_SNPs.normed_vcf,
        bcftools_vcf_tbi = Norm_SNPs.normed_vcf_tbi,
        prefix = prefix
    }

    ##### add merge small + sv vcf step #######
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamplesSVs { input:
        vcf_input = HP_SV.phased_sv_vcf,
        tbi_input = HP_SV.phased_sv_vcf_tbi,
        pref = prefix
    }
    call ConcatVCFs { input:
        bcftools_small_vcf = Filter_SNPs.filtered_vcf,
        bcftools_small_vcf_tbi = Filter_SNPs.filtered_vcf_tbi,
        bcftools_sv_vcf = MergeAcrossSamplesSVs.merged_vcf,
        bcftools_sv_vcf_tbi = MergeAcrossSamplesSVs.merged_tbi,
        prefix = prefix
    }


    ################################
    call StatPhase.Shapeit4 as Shapeit4scaffold { input:
        vcf_input = ConcatVCFs.concat_vcf,
        vcf_index = ConcatVCFs.concat_vcf_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = region,
        num_threads = num_t
    }
    call FF.FinalizeToFile as Finalizescaffold {
        input: outdir = gcs_out_root_dir, file = Shapeit4scaffold.scaffold_vcf
    }
    ##### phase structural variants

    # call StatPhase.Shapeit4_phaseSVs as Shapeit4SVphase { input:
    #     vcf_input = MergeAcrossSamplesSVs.merged_vcf,
    #     vcf_index = MergeAcrossSamplesSVs.merged_tbi,
    #     scaffold_vcf = Shapeit4scaffold.scaffold_vcf,
    #     mappingfile = genetic_mapping_dict[chromosome],
    #     region = region,
    #     num_threads = num_t
    # }
    # call FF.FinalizeToFile as FinalizeSVs {
    #     input: outdir = gcs_out_root_dir, file = Shapeit4SVphase.final_phased_vcf
    # }

    output{
        File snp_shapeit4scaffold = Shapeit4scaffold.scaffold_vcf
        # File sv_shapeit4scaffold = Shapeit4SVphase.final_phased_vcf
        
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


task SplitMultiallelicCalls {

    meta {
        description: "using bcftools to splic multiallelic calls"
    }

    parameter_meta {
        bcftools_merged_vcf: "bcftools merged vcf files"
    }

    input {
        File bcftools_merged_vcf
        File bcftools_merged_vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail
        bcftools norm -m -any ~{bcftools_merged_vcf} | bgzip > ~{prefix}.normed.vcf.gz
        tabix -p vcf ~{prefix}.normed.vcf.gz
        
    >>>

    output {
        File normed_vcf = " ~{prefix}.normed.vcf.gz"
        File normed_vcf_tbi = "~{prefix}.normed.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }

}

task FilterSmallVariants {

    meta {
        description: "using bcftools to filter small variants by F-missing and MAC"
    }

    parameter_meta {
        bcftools_vcf: "bcftools merged vcf files"
    }

    input {
        File bcftools_vcf
        File bcftools_vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail
        bcftools view -i 'F_MISSING < 0.05 & MAC>=2' ~{bcftools_vcf} | bgzip > ~{prefix}.filtered.vcf.gz
        tabix -p vcf ~{prefix}.filtered.vcf.gz
        
    >>>

    output {
        File filtered_vcf = "~{prefix}.filtered.vcf.gz"
        File filtered_vcf_tbi = "~{prefix}.filtered.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }

}

task ConcatVCFs {

    meta {
        description: "using bcftools to concatenate small variants and SVs"
    }

    parameter_meta {
        bcftools_small_vcf: "bcftools merged vcf files"
    }

    input {
        File bcftools_small_vcf
        File bcftools_small_vcf_tbi
        File bcftools_sv_vcf
        File bcftools_sv_vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail
        bcftools concat ~{bcftools_small_vcf} ~{bcftools_sv_vcf} -Oz -o ~{prefix}_integrated.vcf.gz
        bcftools sort ~{prefix}_integrated.vcf.gz -O z -o ~{prefix}_integrated_sorted.vcf.gz
        tabix -p vcf ~{prefix}_integrated_sorted.vcf.gz
        
    >>>

    output {
        File concat_vcf = "~{prefix}_integrated_sorted.vcf.gz"
        File concat_vcf_tbi = "~{prefix}_integrated_sorted.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }

}
