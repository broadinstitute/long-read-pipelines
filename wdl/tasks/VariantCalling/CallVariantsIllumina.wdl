version 1.0

import "../Utility/Utils.wdl"
import "../Utility/VariantUtils.wdl"
import "../VariantCalling/DeepVariant.wdl" as DV

workflow CallVariants {
    meta {
        description: "A workflow for calling small variants from an Illumina BAM file."
    }
    input {
        File bam
        File bai

        String prefix
        String sample_id

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        Boolean call_small_variants
        Boolean call_small_vars_on_mitochondria = true

        Boolean run_dv_pepper_analysis
        Int? dvp_threads
        Int? dvp_memory

        String mito_contig = "chrM"
        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }

    ######################################################################
    # Block for small variants handling
    ######################################################################

    call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}

    # todo: merge the two scattering scheme into a better one
    if (call_small_variants) {
        # Scatter by chromosome
        Array[String] use_filter = if (call_small_vars_on_mitochondria) then contigs_names_to_ignore else flatten([[mito_contig], contigs_names_to_ignore])
        call Utils.MakeChrIntervalList as SmallVariantsScatterPrepp {
            input:
                ref_dict = ref_dict,
                filter = use_filter
        }

        # size-balanced scatter
        if (run_dv_pepper_analysis) {
            call DV.DeepVariant {
                input:
                    bam           = bam,
                    bai           = bai,
                    ref_fasta     = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,

                    pepper_threads = select_first([dvp_threads]),
                    pepper_memory  = select_first([dvp_memory]),
                    dv_threads = select_first([dvp_threads]),
                    dv_memory  = select_first([dvp_memory]),
                    zones = arbitrary.zone_string
            }
        }
    }

    output {
        File? dvp_g_vcf = DeepVariant.gVCF
        File? dvp_g_tbi = DeepVariant.gVCF_tbi
        File? dvp_vcf = DeepVariant.VCF
        File? dvp_tbi = DeepVariant.VCF_tbi
    }
}
