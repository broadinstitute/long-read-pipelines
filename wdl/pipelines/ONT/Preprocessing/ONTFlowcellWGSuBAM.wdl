version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/QC/FPCheckAoU.wdl" as FPCheck
import "../../../tasks/Visualization/NanoPlot.wdl" as NP
import "../../../tasks/Alignment/AlignONTWGSuBAM.wdl" as Alignment

workflow ONTFlowcellWGSuBAM {

    meta {
        description: "Align one ONT flowcell's reads to a reference genome. Note this is NOT suitable for flowcells that have multiple basecall directories."
    }

    parameter_meta {
        # input
        uBAM:      "GCS path to the unaligned BAM of the flowcell. Note that it's assumed the infamous duplicated read issue has been handled."
        flowcell:  "flowcell ID"
        bam_SM_ID: "sample name to use in the aligned BAM output."
        uBAM_tags_to_preserve: "SAM tags to transfer to the alignment bam, from the input BAM. Please be careful with this: methylation tags MM, ML are usually what you want to keep at minimum."

        turn_off_fingperprint_check: "Please turn off fingerprint check if the reference is not GRCh38."
        fingerprint_store:           "Bucket name and prefix (gs://...) storing the fingerprinting VCFs"
        sample_id_at_store: "Name of the sample at the fingerprint store."

        aln_disk_type: "An optimization specifying which type of disk to use for the minimap2 alignment step."

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        # outputs
        alignment_metrics_tar_gz : "A tar.gz file holding the custom alignment metrics."

        fp_status : "A summary string on the fingerprint check result, value is one of [PASS, FAIL, BORDERLINE]."
        fp_lod_expected_sample : "An LOD score assuming the BAM is the same sample as the FP VCF, i.e. BAM sourced from the 'expected sample'."
        fingerprint_detail_tar_gz : "A tar.gz file holding the fingerprinting details."
    }

    input {
        File uBAM
        String flowcell
        String bam_SM_ID
        Array[String] uBAM_tags_to_preserve = ['MM', 'ML']

        String aln_disk_type = 'SSD'

        String? fingerprint_store
        String? sample_id_at_store
        Boolean turn_off_fingperprint_check = false

        File ref_map_file

        String gcs_out_root_dir
    }

    output {
        # Aligned BAM file
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        Float wgs_cov = ALN.wgs_cov

        # Input unaligned read metrics
        Map[String, Float] uBAM_metrics = NanoPlotFromUBam.stats_map
        File uBAM_metrics_tar_gz = FinalizeInputMetrics.gcs_path

        # Aligned read stats
        Map[String, Float] aln_metrics = ALN.aln_summary
        File alignment_metrics_tar_gz = FinalizeAlnMetrics.gcs_path

        String? fingerprint_check_results = fingerprint_check_status
        Float? fingerprint_check_LOD = FPCheckAoU.lod_expected_sample
        File? fingerprint_check_tar_gz = FinalizeFPDetails.gcs_path

        String last_processing_date = today.yyyy_mm_dd
    }

    ###################################################################################
    # arg validation
    if (!turn_off_fingperprint_check) {
        Boolean fp_info_defined = defined(fingerprint_store) && defined(sample_id_at_store)
        if (!fp_info_defined) {
            call Utils.StopWorkflow as FPInputInvalid { input: reason = "Must provide fingerprint_store and sample_id_at_store when performing fingerprint checks."}
        }
    }

    ###################################################################################
    # prep work
    String workflow_name = "ONTFlowcellWGSuBAM"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{flowcell}"
    String outdir_aln     = outdir + '/alignments'
    String outdir_metrics = outdir + '/metrics'

    Float skip_fingerprint_cov_threshold = 0.5  # if coverage is below this threshold, don't bother with fingerprint check

    ###################################################################################
    # some input metric collection
    call NP.NanoPlotFromUBam { input: uBAM = uBAM }
    call GU.TarGZFiles as PackUbamMetrics {
        input: files = NanoPlotFromUBam.plots, name = 'uBAM.metrics'
    }
    call FF.FinalizeToFile as FinalizeInputMetrics { input: outdir = outdir_metrics, file = PackUbamMetrics.you_got_it }
    ###################################################################################
    # align & save
    call Alignment.AlignONTWGSuBAM as ALN {
        input:
            uBAM = uBAM,
            uBAM_tags_to_preserve = uBAM_tags_to_preserve,
            bam_sample_name = bam_SM_ID,
            flowcell = flowcell,
            ref_map_file = ref_map_file,
            aln_disk_type = aln_disk_type
    }
    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = outdir_aln, file = ALN.aligned_bam }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = outdir_aln, file = ALN.aligned_bai }
    call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = outdir_metrics, file = ALN.alignment_metrics_tar_gz }

    ###################################################################################
    # fingerprint & save
    if (!turn_off_fingperprint_check) {
        if (ALN.wgs_cov > skip_fingerprint_cov_threshold) { # only check when requested && coverage isn't hopeless
            Map[String, String] ref_map = read_map(ref_map_file)
            call FPCheck.FPCheckAoU {
                input:
                    aligned_bam = ALN.aligned_bam,
                    aligned_bai = ALN.aligned_bai,
                    fp_store = select_first([fingerprint_store]),
                    sample_id_at_store = select_first([sample_id_at_store]),
                    ref_specific_haplotype_map = ref_map['haplotype_map']
            }
            call GU.TarGZFiles as PackFPResult {
                input: files = [FPCheckAoU.fingerprint_summary, FPCheckAoU.fingerprint_details], name = 'fingerprint_check.summary_and_details'
            }
            call FF.FinalizeToFile as FinalizeFPDetails { input: outdir = outdir_metrics, file = PackFPResult.you_got_it }
        }

        # FAIL the bam if coverage is below a certain threshold (not much useful data anyway)
        String fingerprint_check_status = select_first([FPCheckAoU.FP_status, "FAIL"])
    }

    ###################################################################################
    call GU.GetTodayDate as today {}
}
