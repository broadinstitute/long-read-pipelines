version 1.0

import "../../../tasks/Utility/SRUtils.wdl" as SRUTIL
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/QC/FastQC.wdl" as FastQC
import "../../TechAgnostic/Utility/RemoveSingleOrganismContamination.wdl" as DECONTAMINATE
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/QC/QCAssessment.wdl" as QCAssessment
workflow SRFlowcell_Simplified {

    meta {
        author: "Jonn Smith"
        description: "This workflow preprocesses short read flowcell data in preparation for variant calling.  This workflow contains the following steps: 1) Sam -> Fastq (if necessary), 2) Alignment to reference with bwa-mem2 (https://github.com/bwa-mem2/bwa-mem2), 3) Mark Duplicate reads, 4) Recalibrate base quality scores."
    }

    parameter_meta {
        bam: "Bam file containing reads to align and process.  `bai` must be defined if this argument is.  This argument and `bai` are mutually exclusive with `fq_end1` and `fq_end2`"
        bai: "Index for `bam`.  `bam` must be defined if this argument is.  This argument and `bam` are mutually exclusive with `fq_end1` and `fq_end2`"

        fq_end1: "FASTQ file containing end 1 of the short read data to process.  `fq_end2` must be defined if this argument is.  This argument and `fq_end2` are mutually exclusive with `bam` and `bai`"
        fq_end2: "FASTQ file containing end 2 of the short read data to process.  `fq_end1` must be defined if this argument is.  This argument and `fq_end1` are mutually exclusive with `bam` and `bai`"

        SM: "Sample name for the given bam file."
        LM: "Library name for the given bam file."

        ref_map_file:  "Reference map file for the primary reference sequence and auxillary file locations."
        contaminant_ref_name:  "Name for the contaminant reference."
        contaminant_ref_map_file:  "Reference map file for the contaminant reference sequence and auxillary file locations."

        dir_prefix: "Directory prefix to use for finalized location."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."

        coverage_bed_file: "BED file containing regions to calculate coverage metrics for."

        perform_BQSR: "If true, will perform Base Quality Score Recalibration.  If false will not recalibrate base qualities."

        perform_mark_duplicates: "If true, will run MarkDuplicates on the raw reads from the sequencer."

        DEBUG_MODE: "If true, will add extra logging and extra debugging outputs."

        platform: "Platform on which the sample for the given bam file was sequenced."
    }

    input {
        File? bam
        File? bai

        File? fq_end1
        File? fq_end2

        String SM
        String LB

        File ref_map_file
        String? contaminant_ref_name
        File? contaminant_ref_map_file

        String dir_prefix

        String? gcs_out_root_dir

        File? coverage_bed_file

        Boolean perform_BQSR = true

        Boolean perform_mark_duplicates = true

        Boolean DEBUG_MODE = false

        String platform = "illumina"
    }

    ####################################
    #     _____         _
    #    |_   _|_ _ ___| | _____
    #      | |/ _` / __| |/ / __|
    #      | | (_| \__ \   <\__ \
    #      |_|\__,_|___/_|\_\___/
    #
    ####################################

    # Get ref info:
    Map[String, String] ref_map = read_map(ref_map_file)

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    if (defined(bam)) {
        # Convert the given bam to a uBAM (needed for previous aligned data):
        call SRUTIL.RevertSam as t_002_RevertSam {
            input:
                input_bam = select_first([bam]),
                prefix = SM + ".revertSam"
        }

        # Convert input SAM/BAM to FASTQ:
        call SRUTIL.BamToFq as t_003_Bam2Fastq {
            input:
                bam = t_002_RevertSam.bam,
                prefix = SM
        }

        call Utils.GetRawReadGroup as t_004_GetRawReadGroup { input: gcs_bam_path = select_first([bam]) }
    }

    # Perform decontamination if we have the contaminant reference:
    if (defined(contaminant_ref_map_file)) {

        # Call our sub-workflow for decontamination:
        # NOTE: We don't need to be too concerned with the finalization info.
        #       This will be partially filled in by the WDL itself, so we can pass the same inputs for
        #       these things here (e.g. `dir_prefix`):
        call DECONTAMINATE.RemoveSingleOrganismContamination as t_005_DecontaminateSample {
            input:
                fq_end1 = select_first([fq_end1, t_003_Bam2Fastq.fq_end1]),
                fq_end2 = select_first([fq_end2, t_003_Bam2Fastq.fq_end2]),

                SM = SM,
                LB = LB,
                platform = platform,

                contaminant_ref_name = select_first([contaminant_ref_name]),
                contaminant_ref_map_file = select_first([contaminant_ref_map_file]),

                dir_prefix = dir_prefix,
                gcs_out_root_dir = gcs_out_root_dir
        }
    }

    File fq_e1 = select_first([t_005_DecontaminateSample.decontaminated_fq1, fq_end1, t_003_Bam2Fastq.fq_end1])
    File fq_e2 = select_first([t_005_DecontaminateSample.decontaminated_fq2, fq_end2, t_003_Bam2Fastq.fq_end2])

    String RG = select_first([t_004_GetRawReadGroup.rg, "@RG\tID:" + SM + "_" + LB + "\tPL:" + platform + "\tLB:" + LB + "\tSM:" + SM])

    # Align reads to reference with BWA-MEM2:
    call SRUTIL.BwaMem2 as t_006_AlignReads {
        input:
            fq_end1 = fq_e1,
            fq_end2 = fq_e2,
            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],
            ref_0123 = ref_map["0123"],
            ref_amb = ref_map["amb"],
            ref_ann = ref_map["ann"],
            ref_bwt = ref_map["bwt"],
            ref_pac = ref_map["pac"],
            mark_short_splits_as_secondary = true,
            read_group = RG,
            prefix = SM + ".aligned"
    }

    if (defined(bam)) {
        # Merge aligned reads and unaligned reads:
        call SRUTIL.MergeBamAlignment as t_007_MergeBamAlignment {
            input:
                aligned_bam = t_006_AlignReads.bam,
                unaligned_bam = select_first([t_002_RevertSam.bam]),
                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],
                prefix = SM + ".aligned.merged"
        }
    }

    File merged_bam = select_first([t_007_MergeBamAlignment.bam, t_006_AlignReads.bam])

    if (perform_mark_duplicates) {
        # Mark Duplicates
        call SRUTIL.MarkDuplicatesAndSort as t_008_MarkDuplicates {
            input:
                input_bam = merged_bam,
                prefix = SM + ".aligned.merged.markDuplicates.sorted"
        }
    }
    if (!perform_mark_duplicates) {
        # Sort merged bam:
        call Utils.SortSam as t_009_SortBam {
            input:
                input_bam = merged_bam,
                prefix = basename(merged_bam, ".bam") + ".sorted"
        }
    }
    File pre_bqsr_bam = select_first([t_008_MarkDuplicates.bam, t_009_SortBam.output_bam])
    File pre_bqsr_bam_index = select_first([t_008_MarkDuplicates.bai, t_009_SortBam.output_bam_index])

#    TODO: Add Fingerprinting?

    if (perform_BQSR) {
        # Recalibrate Base Scores:
        call SRUTIL.RunBaseRecalibratorAndApplyBQSR as t_010_RunBaseRecalibratorAndApplyBQSR {
            input:
                input_bam = pre_bqsr_bam,
                input_bam_index = pre_bqsr_bam_index,

                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],

                known_sites_vcf = ref_map["known_sites_vcf"],
                known_sites_index = ref_map["known_sites_index"],

                prefix = SM + ".aligned.merged.markDuplicates.sorted.BQSR"
        }
    }

    File final_bam = select_first([t_010_RunBaseRecalibratorAndApplyBQSR.recalibrated_bam, pre_bqsr_bam])
    File final_bai = select_first([t_010_RunBaseRecalibratorAndApplyBQSR.recalibrated_bai, pre_bqsr_bam_index])

    #############################################
    #      __  __      _        _
    #     |  \/  | ___| |_ _ __(_) ___ ___
    #     | |\/| |/ _ \ __| '__| |/ __/ __|
    #     | |  | |  __/ |_| |  | | (__\__ \
    #     |_|  |_|\___|\__|_|  |_|\___|___/
    #
    #############################################

    call AM.SamStatsMap as t_012_SamStats {
        input:
            bam = final_bam
    }

    call FastQC.FastQC as t_013_FastQC { input: bam = final_bam, bai = final_bai }
    call Utils.ComputeGenomeLength as t_014_ComputeGenomeLength { input: fasta = ref_map['fasta'] }
    call SRUTIL.ComputeBamStats as t_015_ComputeBamStats { input: bam_file = final_bam }

    # Collect stats on aligned reads:
    call SRUTIL.ComputeBamStats as t_016_ComputeBamStatsQ5 { input: bam_file  = final_bam, qual_threshold = 5 }
    call SRUTIL.ComputeBamStats as t_017_ComputeBamStatsQ10 { input: bam_file = final_bam, qual_threshold = 10 }
    call SRUTIL.ComputeBamStats as t_018_ComputeBamStatsQ20 { input: bam_file = final_bam, qual_threshold = 20 }
    call SRUTIL.ComputeBamStats as t_019_ComputeBamStatsQ30 { input: bam_file = final_bam, qual_threshold = 30 }
    call SRUTIL.ComputeBamStats as t_020_ComputeBamStatsQ40 { input: bam_file = final_bam, qual_threshold = 40 }

    call AM.AlignedMetrics as t_021_PerFlowcellMetrics {
        input:
            aligned_bam    = final_bam,
            aligned_bai    = final_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = metrics_dir
    }

    call AM.CallableLoci as t_022_CallableLoci {
        input:
            bam_file = final_bam,
            bam_index = final_bai,
            ref_fasta = ref_map['fasta'],
            ref_fasta_index = ref_map['fai'],
            ref_dict = ref_map['dict'],
            min_depth = 5,
            prefix = SM
    }

    #############################################
    #       ___   ____   ____                    __  _____     _ _ 
    #     / _ \ / ___| |  _ \ __ _ ___ ___     / / |  ___|_ _(_) |
    #    | | | | |     | |_) / _` / __/ __|   / /  | |_ / _` | | |
    #    | |_| | |___  |  __/ (_| \__ \__ \  / /   |  _| (_| | | |
    #     \__\_\\____| |_|   \__,_|___/___/ /_/    |_|  \__,_|_|_|
    #    
    #############################################             
    # Only create QC pass/fail metrics if we have a coverage bed file:
    if (defined(coverage_bed_file)) {
        call AM.MosDepthOverBed as t_023_MosDepthOverBed {
            input:
                bam = final_bam,
                bai = final_bai,
                bed = select_first([coverage_bed_file])
        }

        call QCAssessment.AssessQualityMetrics as t_024_AssessQualityMetrics {
            input:
                callable_loci_summary = t_022_CallableLoci.callable_loci_summary,
                mosdepth_region_bed = t_023_MosDepthOverBed.regions,
                prefix = SM
        }
    }

    ############################################
    #      _____ _             _ _
    #     |  ___(_)_ __   __ _| (_)_______
    #     | |_  | | '_ \ / _` | | |_  / _ \
    #     |  _| | | | | | (_| | | |/ /  __/
    #     |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ############################################

    if (defined(gcs_out_root_dir)) {
        # Create an outdir:
        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])

        String outdir = if DEBUG_MODE then sub(concrete_gcs_out_root_dir, "/$", "") + "/SRFlowcell_Simplified/~{dir_prefix}/" + t_001_WdlExecutionStartTimestamp.timestamp_string else sub(concrete_gcs_out_root_dir, "/$", "") + "/SRFlowcell/~{dir_prefix}"
        String reads_dir = outdir + "/reads"
        String unaligned_reads_dir = outdir + "/reads/unaligned"
        String aligned_reads_dir = outdir + "/reads/aligned"
        String metrics_dir = outdir + "/metrics"

        File keyfile = t_022_CallableLoci.callable_loci_summary

        # Finalize our unaligned reads first:
        call FF.FinalizeToDir as t_022_FinalizeUnalignedFastqReads {
            input:
                outdir = unaligned_reads_dir,
                files =
                [
                    fq_e1,
                    fq_e2,
                ],
                keyfile = keyfile
        }
        if (defined(bam)) {
            call FF.FinalizeToDir as t_023_FinalizeUnalignedReadsFromBam {
                input:
                    outdir = unaligned_reads_dir,
                    files = select_all(
                    [
                        bam,
                        bai,
                        t_003_Bam2Fastq.fq_unpaired,
                    ]),
                    keyfile = keyfile
            }
        }

        call FF.FinalizeToDir as t_024_FinalizeAlignedReads {
            input:
                outdir = aligned_reads_dir,
                files =
                [
                    t_006_AlignReads.bam,
                    merged_bam,
                    pre_bqsr_bam,
                    pre_bqsr_bam_index,
                ],
                keyfile = keyfile
        }

        call FF.FinalizeToFile as t_025_FinalizeAlignedBam {
            input:
                outdir = aligned_reads_dir,
                file = final_bam,
                keyfile = keyfile
        }

        call FF.FinalizeToFile as t_026_FinalizeAlignedBai {
            input:
                outdir = aligned_reads_dir,
                file = final_bai,
                keyfile = keyfile
        }

        # Finalize our metrics:
        call FF.FinalizeToDir as t_027_FinalizeMetrics {
            input:
                outdir = metrics_dir,
                files =
                [
                    t_012_SamStats.sam_stats,
                    t_015_ComputeBamStats.results_file,
                    t_016_ComputeBamStatsQ5.results_file,
                    t_017_ComputeBamStatsQ10.results_file,
                    t_018_ComputeBamStatsQ20.results_file,
                    t_019_ComputeBamStatsQ30.results_file,
                    t_020_ComputeBamStatsQ40.results_file,
                ],
                keyfile = keyfile
        }

        # Finalize BQSR Metrics if it was run:
        if (perform_mark_duplicates) {
            call FF.FinalizeToDir as t_028_FinalizeMarkDuplicatesMetrics {
                input:
                    outdir = metrics_dir,
                    files = select_all([t_008_MarkDuplicates.metrics]),
                    keyfile = keyfile
            }
        }

        # Finalize BQSR Metrics if it was run:
        if (perform_BQSR) {
            call FF.FinalizeToDir as t_029_FinalizeBQSRMetrics {
                input:
                    outdir = metrics_dir,
                    files = select_all([t_010_RunBaseRecalibratorAndApplyBQSR.recalibration_report]),
                    keyfile = keyfile
            }
        }

        call FF.FinalizeToFile as t_030_FinalizeFastQCReport {
            input:
                outdir = metrics_dir,
                file = t_013_FastQC.report
        }

        # Prep a few files for output:
        File fq1_o = unaligned_reads_dir + "/" + basename(fq_e1)
        File fq2_o = unaligned_reads_dir + "/" + basename(fq_e2)
        if (defined(bam)) {
            File unaligned_bam_o = unaligned_reads_dir + "/" + basename(select_first([bam]))
            File unaligned_bai_o = unaligned_reads_dir + "/" + basename(select_first([bai]))
            File fqboup = unaligned_reads_dir + "/" + basename(select_first([t_005_DecontaminateSample.decontaminated_unpaired, t_003_Bam2Fastq.fq_unpaired]))
        }

        call FF.FinalizeToFile as t_031_FinalizeCallableLociSummary {
            input:
                outdir = metrics_dir,
                file = t_022_CallableLoci.callable_loci_summary,
                keyfile = keyfile
        }

        call FF.FinalizeToFile as t_032_FinalizeCallableLociBed {
            input:
                outdir = metrics_dir,
                file = t_022_CallableLoci.callable_loci_bed,
                keyfile = keyfile
        }
    }

    # Prep some output values before the output block:
    Float raw_est_fold_cov_value = if t_014_ComputeGenomeLength.length != 0 then t_015_ComputeBamStats.results['bases']/t_014_ComputeGenomeLength.length else 0.0
    Float aligned_frac_bases_value = if t_012_SamStats.stats_map['total_length'] != 0 then t_012_SamStats.stats_map['bases_mapped']/t_012_SamStats.stats_map['total_length'] else 0.0
    Float aligned_est_fold_cov_value = if t_014_ComputeGenomeLength.length != 0 then t_012_SamStats.stats_map['bases_mapped']/t_014_ComputeGenomeLength.length else 0.0
    Float average_identity_value = if t_012_SamStats.stats_map['bases_mapped'] != 0 then 100.0 - (100.0*t_012_SamStats.stats_map['mismatches']/t_012_SamStats.stats_map['bases_mapped']) else 0.0

    ############################################
    #      ___        _               _
    #     / _ \ _   _| |_ _ __  _   _| |_
    #    | | | | | | | __| '_ \| | | | __|
    #    | |_| | |_| | |_| |_) | |_| | |_
    #     \___/ \__,_|\__| .__/ \__,_|\__|
    #                    |_|
    ############################################
    output {
        # Unaligned reads
        File fq1 = select_first([fq1_o, fq_e1])
        File fq2 = select_first([fq2_o, fq_e1])
        File? fq_unpaired = fqboup

        # Unaligned BAM file
        File? unaligned_bam = unaligned_bam_o
        File? unaligned_bai = unaligned_bai_o

        # Aligned BAM file
        File aligned_bam = select_first([t_025_FinalizeAlignedBam.gcs_path, final_bam])
        File aligned_bai = select_first([t_026_FinalizeAlignedBai.gcs_path, final_bai])

        # Contaminated BAM file:
        # TODO: This will need to be fixed for optional finalization:
        File? contaminated_bam = t_005_DecontaminateSample.contaminated_bam
        File? contaminated_bam_index = t_005_DecontaminateSample.contaminated_bam_index

        ##############################
        # Statistics:

        # Unaligned read stats
        Float num_reads = t_015_ComputeBamStats.results['reads']
        Float num_bases = t_015_ComputeBamStats.results['bases']
        Float raw_est_fold_cov = raw_est_fold_cov_value

        Float read_length = t_015_ComputeBamStats.results['read_mean']

        Float read_qual_mean = t_015_ComputeBamStats.results['mean_qual']
        Float read_qual_median = t_015_ComputeBamStats.results['median_qual']

        Float num_reads_Q5 = t_016_ComputeBamStatsQ5.results['reads']
        Float num_reads_Q10 = t_017_ComputeBamStatsQ10.results['reads']
        Float num_reads_Q20 = t_018_ComputeBamStatsQ20.results['reads']
        Float num_reads_Q30 = t_019_ComputeBamStatsQ30.results['reads']
        Float num_reads_Q40 = t_020_ComputeBamStatsQ40.results['reads']

        # Aligned read stats
        Float aligned_num_reads = t_013_FastQC.stats_map['number_of_reads']
        Float aligned_num_bases = t_012_SamStats.stats_map['bases_mapped']
        Float aligned_frac_bases = aligned_frac_bases_value
        Float aligned_est_fold_cov = aligned_est_fold_cov_value

        Float aligned_read_length = t_013_FastQC.stats_map['read_length']

        Float insert_size_average = t_012_SamStats.stats_map['insert_size_average']
        Float insert_size_standard_deviation = t_012_SamStats.stats_map['insert_size_standard_deviation']
        Float pct_properly_paired_reads = t_012_SamStats.stats_map['percentage_of_properly_paired_reads_%']

        Float average_identity = average_identity_value

        File fastqc_report = select_first([t_030_FinalizeFastQCReport.gcs_path, t_013_FastQC.report])

        File callable_loci_summary = select_first([t_031_FinalizeCallableLociSummary.gcs_path, t_022_CallableLoci.callable_loci_summary])
        File callable_loci_bed = select_first([t_032_FinalizeCallableLociBed.gcs_path, t_022_CallableLoci.callable_loci_bed])

        String? qc_status = t_024_AssessQualityMetrics.qc_status
        String? qc_message = t_024_AssessQualityMetrics.qc_message
    }
}
