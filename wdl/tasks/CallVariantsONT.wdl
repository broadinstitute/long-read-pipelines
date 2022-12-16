version 1.0

import "Utils.wdl"
import "VariantUtils.wdl"
import "ONTPepper.wdl"

workflow CallVariants {
    meta {
        description: "A workflow for calling small and/or structural variants from an aligned ONT BAM file."
    }
    input {
        File bam
        File bai
        String prefix

        File ref_fasta
        File ref_fasta_fai

        String prefix

        Boolean call_small_variants

        Boolean run_dv_pepper_analysis
        Int? dvp_threads
        Int? dvp_memory
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids
    }

    parameter_meta {

        run_dv_pepper_analysis:  "to turn on DV-Pepper analysis or not (non-trivial increase in cost and runtime)"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; needed only when running DV-Pepper"
        ref_scatter_interval_list_ids:     "A file that gives short IDs to the interval_list files; needed only when running DV-Pepper"
    }

    ######################################################################
    # Block for small variants handling
    ######################################################################

    call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}

    # todo: merge the two scattering scheme into a better one
    if (call_small_variants) {
        # size-balanced scatter
        if (run_dv_pepper_analysis) {
            File scatter_interval_list_ids = select_first([ref_scatter_interval_list_ids])
            File scatter_interval_list_loc = select_first([ref_scatter_interval_list_locator])
            Array[String] interval_list_ids   = read_lines(scatter_interval_list_ids)
            Array[String] interval_list_files = read_lines(scatter_interval_list_loc)
            Array[Pair[String, String]] ided_interval_list_files = zip(interval_list_ids, interval_list_files)

            scatter (pair in ided_interval_list_files) {
                call Utils.ResilientSubsetBam as size_balanced_scatter {
                    input:
                        bam = bam,
                        bai = bai,
                        interval_list_file = pair.right,
                        interval_id = pair.left,
                        prefix = basename(bam, ".bam")
                }

                call ONTPepper.Pepper {
                    input:
                        bam           = size_balanced_scatter.subset_bam,
                        bai           = size_balanced_scatter.subset_bai,
                        ref_fasta     = ref_fasta,
                        ref_fasta_fai = ref_fasta_fai,
                        threads       = select_first([dvp_threads]),
                        memory        = select_first([dvp_memory]),
                        zones = arbitrary.zones
                }
            }

            String dvp_prefix = prefix + ".deepvariant_pepper"

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
                input:
                    vcfs     = Pepper.gVCF,
                    prefix   = dvp_prefix + ".g",
                    ref_fasta_fai = ref_fasta_fai
            }

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantPhasedVCFs {
                input:
                    vcfs     = Pepper.phasedVCF,
                    prefix   = dvp_prefix + ".phased",
                    ref_fasta_fai = ref_fasta_fai
            }

            call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
                input:
                    vcfs     = Pepper.VCF,
                    prefix   = dvp_prefix,
                    ref_fasta_fai = ref_fasta_fai
            }

            call Utils.MergeBams {
                input:
                    bams = Pepper.hap_tagged_bam,
                    prefix = prefix +  ".MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged"
            }
        }
    }

  output {
        File? dvp_phased_vcf = MergeDeepVariantPhasedVCFs.vcf
        File? dvp_phased_tbi = MergeDeepVariantPhasedVCFs.tbi
        File? dvp_g_vcf = MergeDeepVariantGVCFs.vcf
        File? dvp_g_tbi = MergeDeepVariantGVCFs.tbi
        File? dvp_vcf = MergeDeepVariantVCFs.vcf
        File? dvp_tbi = MergeDeepVariantVCFs.tbi
    }
}
