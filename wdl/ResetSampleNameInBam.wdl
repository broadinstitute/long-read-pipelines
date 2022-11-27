version 1.0

import "tasks/Utils.wdl"
import "tasks/Finalize.wdl" as FF
import "tasks/utils/CloudUtils.wdl"
import "tasks/utils/BAMutils.wdl" as BU

workflow ResetSampleNameInBam {
    meta {
        description: "For resetting the sample name in a BAM file to the requested value"
    }

    input {
        File  bam
        File? bai
        String new_sample_name

        Boolean skip_single_sample_bam_check = false
    }

    parameter_meta {
        bam: "BAM to operate on, aligned or not, but is assumed to contain @RG lines in its header holding SM information, and the BAM is for a single sample."
        bai: "Accompanying bai index."
        new_sample_name: "New sample name desired (assumed to be a valid value)."
        skip_single_sample_bam_check: "When false, we skip the safety check that the BAM should hold only a single sample name."
    }

    if (! skip_single_sample_bam_check) {
        call BU.CountReadGroupsAndSamples {input: bam = bam, bai = bai}
        if (CountReadGroupsAndSamples.sm_cnt != 1) {
            call Utils.StopWorkflow { input: reason = "Input bam has " + CountReadGroupsAndSamples.sm_cnt + " samples defined in the header."}
        }
    }
    call Utils.ComputeAllowedLocalSSD as guess {input: intended_gb = 1 + ceil(3 * size(bam, "GiB")) }
    call BU.ReSampleBam as ReSample {input: bam = bam, bai = bai, new_sample_name = new_sample_name, out_prefix = "does_not_matter"}

    call CloudUtils.GetBlobFolder {input: blob = bam}
    String outfolder = GetBlobFolder.bucket_and_folder
    String final_out_prefix = basename(bam, '.bam')
    call FF.FinalizeToFile as ff_bam {input: file = ReSample.reheadered_bam, outdir = outfolder, name = final_out_prefix + ".RH.bam"}
    if (defined(ReSample.reheadered_bai)) {
        call FF.FinalizeToFile as ff_bai {input: file = select_first([ReSample.reheadered_bai]), outdir = outfolder, name = final_out_prefix + ".RH.bai"}
    }

    output {
        File  reheadered_bam = ff_bam.gcs_path
        File? reheadered_bai = ff_bai.gcs_path

        # purely Terra trick
        File  original_bam = bam
        File? original_bai = bai
    }
}
