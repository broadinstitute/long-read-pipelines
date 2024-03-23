version 1.0


import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU

import "../../../structs/ReferenceMetadata.wdl"

import "../../../tasks/VariantCalling/DeepVariant.wdl"

workflow FemaleXYDV {
    input  {
        File bam
        File bai

        String model_type

        File ref_bundle_json_file

        File allosomes_interval_list

        String gcs_out_dir
        String prefix
        Int dv_threads
        Int dv_memory
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    output {
        File xy_fixed_gvcf       = FinalizeVcf.gcs_path
        File xy_fixed_gvcf_index = FinalizeTbi.gcs_path
    }

    HumanReferenceBundle ref_bundle = read_json(ref_bundle_json_file)

    ####################################################################################################################################
    # first, subset to X and Y

    call Utils.ResilientSubsetBam as FirstAttempt { input:
        bam = bam,
        bai = bai,
        interval_list_file = allosomes_interval_list,
        interval_id = "chrX-Y",
        prefix = basename(bam, ".bam"),
    }
    if (FirstAttempt.is_samtools_failed) {
        call Utils.ResilientSubsetBam as SecondAttempt { input:
            bam = bam,
            bai = bai,
            interval_list_file = allosomes_interval_list,
            interval_id = "chrX-Y",
            prefix = basename(bam, ".bam"),
        }
        if (SecondAttempt.is_samtools_failed) {
            call BU.SubsetBamToLocusLocal {
                input:
                    bam = bam,
                    bai = bai,
                    interval_list_file = allosomes_interval_list,
                    interval_id = "chrX-Y",
                    prefix = basename(bam, ".bam")
            }
        }
    }
    File allo_bam = select_first([SubsetBamToLocusLocal.subset_bam, SecondAttempt.subset_bam, FirstAttempt.subset_bam])
    File allo_bai = select_first([SubsetBamToLocusLocal.subset_bai, SecondAttempt.subset_bai, FirstAttempt.subset_bai])

    ####################################################################################################################################
    # then, run DeepVariant on the X and Y only

    call DeepVariant.DV { input:
        bam = allo_bam,
        bai = allo_bai,
        ref_fasta = ref_bundle.fasta,
        ref_fasta_fai = ref_bundle.fai,

        model_type = model_type,

        regions = ["chrX", "chrY"],

        threads = dv_threads,
        memory = dv_memory,
        zones = zones
    }
    ####################################################################################################################################
    # save

    String gcs_variants_out_dir = sub(gcs_out_dir, "/$", "") + "/variants/patch/DV_allosomes"
    call FF.FinalizeToFile as FinalizeVcf  { input: outdir = gcs_variants_out_dir, file = DV.gVCF }
    call FF.FinalizeToFile as FinalizeTbi  { input: outdir = gcs_variants_out_dir, file = DV.gVCF_tbi }
    call FF.FinalizeToFile as FinalizeHtml { input: outdir = gcs_variants_out_dir, file = DV.visual_report_html }
}