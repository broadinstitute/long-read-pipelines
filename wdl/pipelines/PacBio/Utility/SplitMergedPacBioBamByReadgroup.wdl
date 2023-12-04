version 1.0

import "../../TechAgnostic/Utility/SplitBamByReadgroup.wdl" as Major

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/PBUtils.wdl"

workflow SplitMergedPacBioBamByReadgroup {
    meta {
        description: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
    }

    input {
        File input_bam
        File? input_bai

        Boolean unmap_bam
        Boolean convert_to_fq = false
        Boolean disable_primrose_check = false

        String gcs_out_root_dir
    }

    output {
        Map[String, String] rgid_2_bam = WORKHORSE.rgid_2_bam
        Map[String, String] rgid_2_PU  = WORKHORSE.rgid_2_PU
        Map[String, String]? rgid_2_ubam_emptyness = WORKHORSE.rgid_2_ubam_emptyness
        Boolean rgid_2_bam_are_aligned = WORKHORSE.rgid_2_bam_are_aligned
        Map[String, String]? rgid_2_fastq = WORKHORSE.rgid_2_fastq

        String last_processing_date = WORKHORSE.last_processing_date
    }

    parameter_meta {
        input_bam: "BAM to be split by read group; doesn't necessarily need to be aligned or sorted."
        input_bai: "(optional) BAI accompanying the BAM"
        platform:  "long reads platform the BAM was generated on; must be one of [PB, ONT]"
        convert_to_fq: "user option to convert to FASTQ (gz) or not"
        gcs_out_root_dir: "place to store the result files"
        disable_primrose_check: "when true, disable a QC check making sure primrose is run on all readgroups; only do this when you know what you're doing"
    }

    ##############################################################################################################################
    # input BAM most basic QC check
    call BU.GatherBamMetadata {
        input: bam = input_bam
    }
    if ('coordinate' != GatherBamMetadata.sort_order) {
        call Utils.StopWorkflow { input: reason = "Input bam isn't coordinate-sorted, but rather sorted by ~{GatherBamMetadata.sort_order}"  }
    }

    # this guarantees that there are no read groups missing primrose runs
    if (!disable_primrose_check) {
        call PBUtils.VerifyPacBioBamHasAppropriatePrimroseRuns as PrimroseCheck { input: bam = input_bam }
        if (0!= length(PrimroseCheck.readgroups_missing_primrose)) {
            call Utils.StopWorkflow as MissingPrimrose { input: reason = "Input BAM file has some of its read groups missing primrose calls."}
        }
    }

    # delegate
    String workflow_name = 'SplitMergedPacBioBamByReadgroup'
    call Major.SplitBamByReadgroup as WORKHORSE {
        input:
            input_bam = input_bam,
            input_bai = input_bai,

            unmap_bam = unmap_bam,
            convert_to_fq = convert_to_fq,

            validate_output_bams = true,

            gcs_out_root_dir = gcs_out_root_dir,
            override_workflow_name = workflow_name,
            debug_mode = false
    }
}
