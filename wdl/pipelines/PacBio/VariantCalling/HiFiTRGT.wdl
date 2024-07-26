version 1.0

import "../../../tasks/VariantCalling/TRGT.wdl" as TRGT

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Utility/Utils.wdl"

workflow runTRGT {

    meta {
        description: "Uses TRGT to size TRs in a bam file."
    }

    input {
        File input_bam
        File input_bam_bai
        String sex

        String? custom_out_prefix  # if not provided, will use sample name encoded in the BAM as prefix for output files

        File ref_fasta
        File ref_fasta_index

        File repeatCatalog
        String catalog_name  # used in naming, so be linux friendly

        String gcs_out_root_dir
    }

    output {
        File trgt_output_vcf     = FinalizeTrgtVCF.gcs_path
        File trgt_output_vcf_idx = FinalizeTrgtTBI.gcs_path
        # File trgt_output_bam = processWithTRGT.trgt_output_bam
    }

    if (!defined(custom_out_prefix)){
        call Utils.InferSampleName { input: bam = input_bam, bai = input_bam_bai }
    }

    Boolean is_female = 'F'==sex

    call TRGT.processWithTRGT as processWithTRGT { input:
        input_bam = input_bam,
        input_bam_bai = input_bam_bai,
        is_female = is_female,
        outprefix = select_first([custom_out_prefix, InferSampleName.sample_name]),
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        repeatCatalog = repeatCatalog,
        catalog_name = catalog_name
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/TRGT/" + catalog_name
    call FF.FinalizeToFile as FinalizeTrgtVCF { input: outdir = outdir, file = processWithTRGT.trgt_output_vcf     }
    call FF.FinalizeToFile as FinalizeTrgtTBI { input: outdir = outdir, file = processWithTRGT.trgt_output_vcf_idx }
}
