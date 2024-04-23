version 1.0

import "../../../tasks/Phasing/Hiphase.wdl"
import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase

workflow ReadManifestFilesHiphase {
    input {
        File snp_vcf_manifest
        File snp_vcf_tbi_manifest
        File sv_vcf_manifest
        File sv_vcf_tbi_manifest
        File bam_manifest
        File bai_manifest
        File reference_fasta
        File reference_index
        #File genetic_mapping_tsv_for_shapeit4

        String chromosome
        
    }

    call ReadFile as snp_vcf_input { input: input_file = snp_vcf_manifest, chr = chromosome, prefix = "snp_vcf" }
    call ReadFile as snp_tbi_input { input: input_file = snp_vcf_tbi_manifest, chr = chromosome, prefix = "snp_tbif" }
    call ReadFile as sv_vcf_input { input: input_file = sv_vcf_manifest, chr = chromosome, prefix = "sv_vcf" }
    call ReadFile as sv_tbi_input { input: input_file = sv_vcf_tbi_manifest, chr = chromosome, prefix = "sv_tbi" }
    call ReadFile as bam_input { input: input_file = bam_manifest, chr = chromosome, prefix = "bam" }
    call ReadFile as bai_input { input: input_file = bai_manifest, chr = chromosome, prefix = "bai" }

    #Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)
    Int data_length = length(snp_vcf_input.sample)
    Array[Int] indexes= range(data_length)

    scatter (i in indexes) {
        String sample = snp_vcf_input.sample[i]
        String snp_vcf = snp_vcf_input.summary[i]
        String snp_tbi = snp_tbi_input.summary[i]
        String sv_vcf = sv_vcf_input.summary[i]
        # String sv_tbi = sv_tbi_input.summary[i]
        String bam = bam_input.summary[i]
        String bai = bai_input.summary[i]

        call convert as cleanvcf{
            input:
                vcf = sv_vcf,
                samplename = sample
        }

        call SplitMultiallelicCalls as Norm_SNPs { input:
            bcftools_merged_vcf = snp_vcf,
            bcftools_merged_vcf_tbi = snp_tbi,
            prefix = sample
        }

        call Hiphase.HiphaseSVs as HP_SV { input:
            bam = bam,
            bai = bai,
            unphased_snp_vcf = Norm_SNPs.normed_vcf,
            unphased_snp_tbi = Norm_SNPs.normed_vcf_tbi,
            unphased_sv_vcf = cleanvcf.subset_vcf,
            unphased_sv_tbi = cleanvcf.subset_tbi,
            ref_fasta = reference_fasta,
            ref_fasta_fai = reference_index,
            samplename = sample
        }
    }
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
        vcf_input = HP_SV.phased_snp_vcf,
        tbi_input = HP_SV.phased_snp_vcf_tbi,
        pref = chromosome
    }

    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamplesSVs { input:
        vcf_input = HP_SV.phased_sv_vcf,
        tbi_input = HP_SV.phased_sv_vcf_tbi,
        pref = chromosome
    }

}

task ReadFile {
    meta {
        desciption:
        "Parse and convert most critical metrics from the QUAST report on the [primary, H1 and H2] assemblies"
    }

    parameter_meta {
        input_file: "input manifest files"
        chr: "chromosome names"
    }

    input {
        File input_file
        String chr
        String prefix
    }

    output {
        Array[String] summary = read_lines("~{prefix}.path.tsv")
        Array[String] sample = read_lines("~{prefix}.sample.tsv")
    }

    command <<<
        set -eux

        python <<CODE
        import pandas as pd

        df = pd.read_csv("~{input_file}", sep = "\t", index_col = 0)
        dv = df["~{chr}"]
        dv = dv.to_dict()


        with open("~{prefix}.path.tsv", 'w') as outf:
            for sample, path in dv.items():
                outf.write(f"{path}\n")

        with open("~{prefix}.sample.tsv", 'w') as outf:
            for sample, path in dv.items():
                outf.write(f"{sample}\n")
                #outf.write(f"{sample}\t{path}\n")
        CODE
    >>>

    runtime {
        disks: "local-disk 10 HDD"
        docker: "hangsuunc/agg:v1"
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
        cpu: 1
        memory:  "4 GiB"
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
        bcftools norm -m -any ~{bcftools_merged_vcf} | bgzip > "~{prefix}.normed.vcf.gz"
        tabix -p vcf "~{prefix}.normed.vcf.gz"
        
    >>>

    output {
        File normed_vcf = "~{prefix}.normed.vcf.gz"
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
