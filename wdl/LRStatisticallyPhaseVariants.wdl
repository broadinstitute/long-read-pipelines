version 1.0

#######################################################################
## A workflow that performs local ancestry inference in long read data.
#######################################################################

import "tasks/IntervalUtils.wdl"
import "tasks/VariantUtils.wdl"
import "tasks/SHAPEIT5.wdl"
import "tasks/Finalize.wdl" as FF

workflow LRStatisticallyPhaseVariants {
    input {
        File gvcf
        File tbi
        File ref_map_file

        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        gvcf:             "GCS path to gVCF file"
        tbi:              "GCS path to gVCF tbi file"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/LRStatisticallyPhaseVariants/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    call IntervalUtils.GetContigNames { input: ref_dict = ref_map['dict'] }

    scatter (contig_name in [ GetContigNames.contig_names[0] ]) {
        call IntervalUtils.GenerateIntervals as GenerateCommonVariantIntervals {
            input:
                ref_dict        = ref_map['dict'],
                selected_contig = contig_name,
                chunk_bp        = 6000000,
                stride_bp       = 2000000,
                buffer_bp       = 0,
        }

        call VariantUtils.SubsetVCF { input: vcf_gz = gvcf, vcf_tbi = tbi, locus = contig_name }

        call VariantUtils.FillTags { input: vcf_gz = SubsetVCF.subset_vcf, tags = [ 'AC', 'AN' ] }

        scatter (interval in [ GenerateCommonVariantIntervals.intervals[0], GenerateCommonVariantIntervals.intervals[1] ]) {
            call SHAPEIT5.PhaseCommonVariants {
                input:
                    input_vcf  = FillTags.filled_bcf,
                    input_tbi  = FillTags.filled_tbi,
                    filter_maf = 0.005,
                    interval   = interval,
            }
        }

        call SHAPEIT5.LigatePhasedCommonVariants {
            input:
                phased_shard_bcfs = PhaseCommonVariants.common_phased_shard_bcf,
                phased_shard_tbis = PhaseCommonVariants.common_phased_shard_tbi,
                prefix = "out",
        }

        call IntervalUtils.GenerateIntervals as InputRegionIntervals {
            input:
                ref_dict        = ref_map['dict'],
                selected_contig = contig_name,
                chunk_bp        = 6000000,
                stride_bp       = 3000000,
                buffer_bp       = 1500000,
        }

        call IntervalUtils.GenerateIntervals as ScaffoldRegionIntervals {
            input:
                ref_dict        = ref_map['dict'],
                selected_contig = contig_name,
                chunk_bp        = 6000000,
                stride_bp       = 3000000,
                buffer_bp       = 0
        }

        scatter (p in zip([ InputRegionIntervals.intervals[0] ], [ ScaffoldRegionIntervals.intervals[0] ])) {
            call SHAPEIT5.PhaseRareVariants {
                input:
                    input_bcf       = SubsetVCF.subset_vcf,
                    scaffold_bcf    = LigatePhasedCommonVariants.scaffold_bcf,
                    input_region    = p.left,
                    scaffold_region = p.right,
            }
        }

        call SHAPEIT5.ConcatenateVariants as ConcatenateContigVariants {
            input:
                shard_bcfs = PhaseRareVariants.rare_phased_shard_bcf,
                shard_tbis = PhaseRareVariants.rare_phased_shard_tbi,
                prefix = contig_name
        }
    }

    call SHAPEIT5.ConcatenateVariants as ConcatenateAllVariants {
        input:
            shard_bcfs = ConcatenateContigVariants.full_bcf,
            shard_tbis = ConcatenateContigVariants.full_tbi,
            prefix = prefix
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeBCF { input: outdir = outdir, file = ConcatenateAllVariants.full_bcf }
    call FF.FinalizeToFile as FinalizeTBI { input: outdir = outdir, file = ConcatenateAllVariants.full_tbi }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File phased_bcf = FinalizeBCF.gcs_path
        File phased_tbi = FinalizeTBI.gcs_path
    }
}
