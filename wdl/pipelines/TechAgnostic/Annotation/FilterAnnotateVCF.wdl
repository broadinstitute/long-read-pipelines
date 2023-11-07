version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Annotation/VariantFiltration.wdl" as GATK

workflow FilterAnnotateVCF {
    meta {
        desciption: "A simple workflow to apply a filter annotation (to the FILTER column) to an VCF"
    }
    parameter_meta {
        # input
        vcf: "vcf to be annotated"
        out_appendix: "appendix to add to the input VCF's basename"

        filter_exp: "argument to be passed to `--filter-expression` of VariantFiltration in GATK"
        filter_name: "argument to be passed to `--filter-name` of VariantFiltration in GATK"

        gatk_docker_tag: "Tag of the GATK GCR image, hosted at us.gcr.io/broad-gatk/gatk. Example: 4.4.0.0"

        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"

        # output
        annotated_vcf: "the annotated VCF"
    }
    input {
        File vcf
        File tbi
        File ref_map_file
        String out_appendix
        String filter_exp
        String filter_name
        String gatk_docker_tag

        String gcs_out_root_dir
    }
    output {
        File annotated_vcf = FinalizeAVcf.gcs_path
        File annotated_vcf_tbi = FinalizeATbi.gcs_path
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    Map[String, String] ref_map = read_map(ref_map_file)
    call GATK.VariantFiltration as VF { input:
        vcf = vcf,
        tbi = tbi,
        ref_fasta = ref_map['fasta'],
        ref_fai = ref_map['fai'],
        ref_dict = ref_map['dict'],
        out_appendix = out_appendix,
        filter_exp = filter_exp,
        filter_name = filter_name,
        gatk_docker_tag = gatk_docker_tag
    }

    call FF.FinalizeToFile as FinalizeAVcf { input: outdir = outdir, file = VF.annotated_vcf }
    call FF.FinalizeToFile as FinalizeATbi { input: outdir = outdir, file = select_first([VF.annotated_tbi]) }
}
