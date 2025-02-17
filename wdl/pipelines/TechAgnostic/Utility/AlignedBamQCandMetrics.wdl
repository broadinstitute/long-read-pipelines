version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/BAMutils.wdl" as BU

import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/Visualization/NanoPlot.wdl" as NP

import "../../../tasks/QC/FPCheckAoU.wdl" as QC0
import "../../TechAgnostic/Utility/LongReadsContaminationEstimation.wdl" as QC1
import "../../TechAgnostic/Utility/SexCheckNaive.wdl" as QC2
import "../../TechAgnostic/Utility/CountTheBeans.wdl" as QC3

import "SaveFilesToDestination.wdl" as SAVE

workflow Work {
    meta {
        desciption:
        "A workflow that unifies standard human WGS aligned BAM QC checks and metrics collection."
    }

    parameter_meta {

        #########
        # inputs
        tech:
        "The technology used to generate this BAM. Currently, the following values are accepted: [ONT, Sequel, Revio]."

        gcs_out_root_dir:
        "output files will be copied over there"

        cov_bed:
        "An optional BED file on which coverage will be collected (a mean value for each interval)"
        cov_bed_descriptor:
        "A short description of the BED provided for targeted coverage estimation; will be used in naming output files."

        fingerprint_vcf_store:
        "A GCS 'folder' holding fingerprint VCF files"
        fingerprint_sample_id:
        "The ID of the sample supposedly this BAM belongs to; note that the fingerprint VCF is assumed to be located at {fingerprint_vcf_store}/{fingerprint_sample_id}*.vcf(.gz?)"

        vbid2_config_json:
        "A config json to for running the VBID2 contamination estimation sub-workflow; if provided, will trigger the VBID2 sub-workflow for cross-(human)individual contamination estimation."

        expected_sex_type:
        "If provided, triggers sex concordance check. Accepted value: [M, F, NA, na]"

        methyl_tag_check_bam_descriptor:
        "If provided, triggers workflow that collects information on reads that miss MM/ML SAM tags; this is meant to be a short description of the purpose of the BAM (e.g. input, a single readgroup, per sample, etc; doesn't need to be single-file specific); used for saving the reads that miss MM/ML tags."

        #########
        # outputs
        wgs_cov:
        "whole genome mean coverage"

        # nanoplot_summ:
        # "Summary on alignment metrics provided by Nanoplot (todo: study the value of this output)"

        sam_flag_stats:
        "SAM flag stats"

        contamination_est:
        "cross-(human)individual contamination estimation by VerifyBAMID2"

        inferred_sex_info:
        "Inferred sex concordance information if expected sex type is provided"

        methyl_tag_simple_stats:
        "Simple stats on the reads with & without SAM methylation tags (MM/ML)."
        save_methyl_uncalled_reads:
        "If to save the reads without MM/ML tags."

        aBAM_metrics_files:
        "A map where keys are summary-names and values are paths to files generated from the various QC/metrics tasks"
    }

    input {
        File bam
        File bai

        String tech

        File?   cov_bed
        String? cov_bed_descriptor

        String? fingerprint_vcf_store
        String? fingerprint_sample_id

        File? vbid2_config_json
        String? expected_sex_type
        String? methyl_tag_check_bam_descriptor
        Boolean? save_methyl_uncalled_reads

        File ref_map_file
        String disk_type

        String output_prefix # String output_prefix = bam_sample_name + "." + flowcell
        String gcs_out_root_dir
    }

    output {
        Float wgs_cov                     = MosDepthWGS.wgs_cov
        # Map[String, Float] nanoplot_summ  = NanoPlotFromBam.stats_map
        Map[String, Float] sam_flag_stats = ParseFlagStatsJson.qc_pass_reads_SAM_flag_stats

        # fingerprint
        Map[String, String]? fingerprint_check = fp_res

        # contam
        Float? contamination_est = VBID2.contamination_est

        # sex concordance
        Map[String, String]? inferred_sex_info = SexConcordance.inferred_sex_info

        # methyl
        Map[String, String]? methyl_tag_simple_stats = NoMissingBeans.methyl_tag_simple_stats

        # file-based outputs all packed into a finalization map
        Map[String, String] aBAM_metrics_files = FF.result
    }

    String workflow_name = "AlignedBamQCandMetrics"
    String metrics_output_dir = sub(gcs_out_root_dir, "/+$","") + "/~{workflow_name}/~{output_prefix}"

    ###################################################################################
    # arg validation and prep
    ###################################################################################
    Map[String, String] ref_map = read_map(ref_map_file)

    if (defined(fingerprint_vcf_store) != defined(fingerprint_sample_id)) {
        call Utils.StopWorkflow as MisingFingerprintArgs { input:
            reason = "fingerprint_vcf_store and fingerprint_sample_id must be specified together or omitted together"
        }
    }

    if (defined(cov_bed) != defined(cov_bed_descriptor)) {
        call Utils.StopWorkflow as MisingCoverageBEDdescriptor { input:
            reason = "cov_bed and cov_bed_descriptor must be specified together or omitted together"
        }
    }
    ###################################################################################
    # ALWAYS ON QC/METRICS
    ###################################################################################
    ################################
    # coverage
    call AM.MosDepthWGS { input:
        bam = bam, bai = bai, disk_type = disk_type,
        bed = cov_bed, bed_descriptor = cov_bed_descriptor
    }
    FinalizationManifestLine a = object
                                 {files_to_save: [MosDepthWGS.summary_txt],
                                  is_singleton_file: true,
                                  destination: metrics_output_dir,
                                  output_attribute_name: "cov_per_chr"}

    ################################
    # SAM flag stats
    call BU.SamtoolsFlagStats  { input: bam = bam, output_format = 'JSON', disk_type = disk_type }
    call BU.ParseFlagStatsJson { input: sam_flag_stats_json = SamtoolsFlagStats.flag_stats }

    ################################
    # nanoplot
    # call NP.NanoPlotFromBam { input: bam = bam, bai = bai, disk_type = disk_type }
    # FinalizationManifestLine b = object
    #                              {files_to_save: flatten([[NanoPlotFromBam.stats], NanoPlotFromBam.plots]),
    #                               is_singleton_file: false,
    #                               destination: metrics_output_dir + "/nanoplot",
    #                               output_attribute_name: "nanoplot"}

    ###################################################################################
    # OPTIONAL QC/METRICS
    ###################################################################################
    ################################
    # (optional) fingerprint
    if (defined(fingerprint_vcf_store)) {
        call QC0.FPCheckAoU as fingerprint {
            input:
                aligned_bam = bam,
                aligned_bai = bai,
                tech = tech,
                fp_vcf_store = select_first([fingerprint_vcf_store]),
                fp_sample_id = select_first([fingerprint_sample_id]),
                ref_specific_haplotype_map = ref_map['haplotype_map']
        }
        Map[String, String] fp_res = {'status': fingerprint.FP_status,
                                      'LOD': fingerprint.lod_expected_sample}
        String dummy = fingerprint.fingerprint_summary # this supposedly File type output may be a null because the input BAM is teeny
        if ("None"!=dummy) {
            FinalizationManifestLine c = object
                                         {files_to_save: [fingerprint.fingerprint_summary, fingerprint.fingerprint_details],
                                          pack_name: "fingerprint.details.tar.gz",
                                          is_singleton_file: false,
                                          destination: metrics_output_dir,
                                          output_attribute_name: "fingerprint_check"}
        }
    }
    ################################
    # (optional) contamination
    if (defined(vbid2_config_json)) {
        VBID2_config vb_conf = read_json(select_first([vbid2_config_json]))
        call QC1.LongReadsContaminationEstimation as VBID2 { input:
            bam=bam,
            bai=bai,
            ref_map_file=ref_map_file,
            tech = tech,
            gt_sites_bed  = vb_conf.genotyping_sites,
            is_hgdp_sites = vb_conf.is_hgdp_sites,
            is_100k_sites = vb_conf.is_100k_sites,
            disable_baq   = vb_conf.disable_baq,
            max_retries = vb_conf.max_retries,
            disk_type = disk_type,
        }
        # no file to save from contam.est.
    }
    ################################
    # (optional) sex concordance
    if (defined(expected_sex_type)) {
        call QC2.SexCheckNaive as SexConcordance { input:
            bam=bam,
            bai=bai,
            expected_sex_type=select_first([expected_sex_type]),
            mosdepth_summary_txt=MosDepthWGS.summary_txt
        }
        # no file to save from sex concordance
    }
    ################################
    # (optional) verify methylation tags aren't missing
    if (defined(methyl_tag_check_bam_descriptor)) {
        call QC3.CountTheBeans as NoMissingBeans { input:
            bam=bam,
            bai=bai,
            bam_descriptor=select_first([methyl_tag_check_bam_descriptor]),
            gcs_out_root_dir= metrics_output_dir,
            save_read_names_only = !select_first([save_methyl_uncalled_reads, true]),
            use_local_ssd=disk_type=='LOCAL'
        }
        Map[String, String] methyl_out = {"missing_methyl_tag_reads":
                                          NoMissingBeans.methyl_tag_simple_stats['files_holding_reads_without_tags']}
    }

    ###################################################################################
    # save results
    ###################################################################################
    call SAVE.SaveFilestoDestination as FF { input:
        # instructions = select_all([a, b, c]),
        instructions = select_all([a, c]),
        already_finalized = select_all([methyl_out]),
        key_file = select_first(select_all([MosDepthWGS.summary_txt,
                                            # NanoPlotFromBam.stats,
                                            fingerprint.fingerprint_summary,]))
    }
}

# if you want to call me as sub-workflow in a standard run
struct AlignedBamQCnMetricsConfig {
    File?   cov_bed                         # An optional BED file on which coverage will be collected (a mean value for each interval)
    String? cov_bed_descriptor              # A short description of the BED provided for targeted coverage estimation; will be used in naming output files.

    String? fingerprint_vcf_store           # A GCS 'folder' holding fingerprint VCF files

    File?   vbid2_config_json               # A config json to for running the VBID2 contamination estimation sub-workflow;
                                            # if provided, will trigger the VBID2 sub-workflow for cross-(human)individual contamination estimation.

    String? methyl_tag_check_bam_descriptor # A one-word description of the purpose of the BAM
                                            # (e.g. 'rg_post_aln' for "readgroup post alignment", "sm_post_merge" for "sample post merging");
                                            # used for saving the reads that miss MM/ML tags;
                                            # doesn't need to be single-file specific
    Boolean? save_methyl_uncalled_reads     # if true, will save the actual reads that miss the methylation SAM tags (otherwise we only save the read names)
}
