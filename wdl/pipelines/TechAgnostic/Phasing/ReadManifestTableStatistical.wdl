version 1.0

import "../../../tasks/Phasing/Hiphase.wdl"
import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase

workflow ReadManifestFilesHiphase {
    input {
        File snp_merged_vcf
        File snp_merged_vcf_tbi
        File sv_merged_vcf
        File sv_merged_vcf_tbi

        File genetic_mapping_tsv_for_shapeit4

        String chromosome
        
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    call StatPhase.Shapeit4 as Shapeit4scaffold { input:
        vcf_input = snp_merged_vcf,
        vcf_index = snp_merged_vcf_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome
    }

    call StatPhase.Shapeit4_phaseSVs as Shapeit4SVphase { input:
        vcf_input = sv_merged_vcf,
        vcf_index = sv_merged_vcf_tbi,
        scaffold_vcf = Shapeit4scaffold.scaffold_vcf,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome
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
        # dv = dv.to_dict()


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
