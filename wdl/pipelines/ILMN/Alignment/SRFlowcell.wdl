version 1.0

import "../../../tasks/Utility/SRUtils.wdl" as SRUTIL
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/QC/FastQC.wdl" as FastQC
import "../../../tasks/Preprocessing/RemoveSingleOrganismContamination.wdl" as DECONTAMINATE
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRFlowcell {

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
        LB: "Library name for the given bam file."

        ref_map_file:  "Reference map file for the primary reference sequence and auxillary file locations."
        contaminant_ref_name:  "Name for the contaminant reference."
        contaminant_ref_map_file:  "Reference map file for the contaminant reference sequence and auxillary file locations."

        dir_prefix: "Directory prefix to use for finalized location."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs."

        perform_BQSR: "If true, will perform Base Quality Score Recalibration.  If false will not recalibrate base qualities."

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

        String gcs_out_root_dir
        
        Boolean perform_BQSR = true

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

    # Create an outdir:
    String outdir = if DEBUG_MODE then sub(gcs_out_root_dir, "/$", "") + "/SRFlowcell/~{dir_prefix}/" + t_001_WdlExecutionStartTimestamp.timestamp_string else sub(gcs_out_root_dir, "/$", "") + "/SRFlowcell/~{dir_prefix}"
    String reads_dir = outdir + "/reads"
    String unaligned_reads_dir = outdir + "/reads/unaligned"
    String aligned_reads_dir = outdir + "/reads/aligned"
    String metrics_dir = outdir + "/metrics"

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

    # OK, this is inefficient, but let's NOW extract our contaminated reads if we have the info.
    # TODO: Move this into the sections above to make it more efficient.  Specifically where we convert bam -> fastq.
    # TODO: Re-enable this section after decontamination is fixed.  The alignment based method with BWA-MEM doesn't work.  Not clear why, but this does seem somewhat inadequate (simplistic alignment-based strategies).
    # Update: Enabled bowtie2-based decontamination of reads
    if (defined(contaminant_ref_map_file)) {
        # Call our sub-workflow for decontamination:
        # NOTE: We don't need to be too concerned with the finalization info.
        #       This will be partially filled in by the WDL itself, so we can pass the same inputs for
        #       these things here (e.g. `dir_prefix`):
        call DECONTAMINATE.RemoveSingleOrganismContamination as DecontaminateSample {
            input:
                fq_end1 = select_first([fq_end1, t_003_Bam2Fastq.fq_end1]),
                fq_end2 = select_first([fq_end2, t_003_Bam2Fastq.fq_end2]),

                SM = SM,
                LB = LB,
                platform = platform,
                aligner = "bowtie2",

                contaminant_ref_name = select_first([contaminant_ref_name]),
                contaminant_ref_map_file = select_first([contaminant_ref_map_file]),

                dir_prefix = dir_prefix,
                gcs_out_root_dir = gcs_out_root_dir
        }
    }

    File fq_e1 = select_first([DecontaminateSample.decontaminated_fq1, fq_end1, t_003_Bam2Fastq.fq_end1])
    File fq_e2 = select_first([DecontaminateSample.decontaminated_fq2, fq_end2, t_003_Bam2Fastq.fq_end2])

    String RG = select_first([t_004_GetRawReadGroup.rg, "@RG\tID:" + SM + "_" + LB + "\tPL:" + platform + "\tLB:" + LB + "\tSM:" + SM])

    # Align reads to reference with BWA-MEM2:
    call SRUTIL.BwaMem2 as t_005_AlignReads {
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
        call SRUTIL.MergeBamAlignment as t_006_MergeBamAlignment {
            input:
                aligned_bam = t_005_AlignReads.bam,
                unaligned_bam = select_first([t_002_RevertSam.bam]),
                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],
                prefix = SM + ".aligned.merged"
        }
    }

    File merged_bam = select_first([t_005_AlignReads.bam, t_006_MergeBamAlignment.bam])

    # Mark Duplicates
    call SRUTIL.MarkDuplicates as t_007_MarkDuplicates {
        input:
            input_bam = merged_bam,
            prefix = SM + ".aligned.merged.markDuplicates"
    }

    # Sort Duplicate Marked Bam:
    call Utils.SortSam as t_008_SortAlignedDuplicateMarkedBam {
        input:
            input_bam = t_007_MarkDuplicates.bam,
            output_bam_basename = SM + ".aligned.merged.markDuplicates.sorted",
            compression_level = 2
    }

#    TODO: Add Fingerprinting?

    if (perform_BQSR) {
        # Recalibrate Base Scores:
        call SRUTIL.BaseRecalibrator as t_009_BaseRecalibrator {
            input:
                input_bam = t_008_SortAlignedDuplicateMarkedBam.output_bam,
                input_bam_index = t_008_SortAlignedDuplicateMarkedBam.output_bam_index,

                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],

                known_sites_vcf = ref_map["known_sites_vcf"],
                known_sites_index = ref_map["known_sites_index"],

                prefix = SM + ".baseRecalibratorReport"
        }

        call SRUTIL.ApplyBQSR as t_010_ApplyBQSR {
            input:
                input_bam = t_008_SortAlignedDuplicateMarkedBam.output_bam,
                input_bam_index = t_008_SortAlignedDuplicateMarkedBam.output_bam,

                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],

                recalibration_report = t_009_BaseRecalibrator.recalibration_report,

                prefix = SM + ".aligned.merged.markDuplicates.sorted.BQSR"
        }
    }

    File final_bam = select_first([t_010_ApplyBQSR.recalibrated_bam, t_008_SortAlignedDuplicateMarkedBam.output_bam])
    File final_bai = select_first([t_010_ApplyBQSR.recalibrated_bai, t_008_SortAlignedDuplicateMarkedBam.output_bam_index])

    #############################################
    #      __  __      _        _
    #     |  \/  | ___| |_ _ __(_) ___ ___
    #     | |\/| |/ _ \ __| '__| |/ __/ __|
    #     | |  | |  __/ |_| |  | | (__\__ \
    #     |_|  |_|\___|\__|_|  |_|\___|___/
    #
    #############################################

    call AM.SamStatsMap as t_011_SamStats {
        input:
            bam = final_bam
    }

    call FastQC.FastQC as t_012_FastQC { input: bam = final_bam, bai = final_bai }
    call Utils.ComputeGenomeLength as t_013_ComputeGenomeLength { input: fasta = ref_map['fasta'] }
    call SRUTIL.ComputeBamStats as t_014_ComputeBamStats { input: bam_file = final_bam }

    # Collect stats on aligned reads:
    call SRUTIL.ComputeBamStats as t_015_ComputeBamStatsQ5 { input: bam_file  = final_bam, qual_threshold = 5 }
    call SRUTIL.ComputeBamStats as t_016_ComputeBamStatsQ7 { input: bam_file  = final_bam, qual_threshold = 7 }
    call SRUTIL.ComputeBamStats as t_017_ComputeBamStatsQ10 { input: bam_file = final_bam, qual_threshold = 10 }
    call SRUTIL.ComputeBamStats as t_018_ComputeBamStatsQ12 { input: bam_file = final_bam, qual_threshold = 12 }
    call SRUTIL.ComputeBamStats as t_019_ComputeBamStatsQ15 { input: bam_file = final_bam, qual_threshold = 15 }

    call AM.AlignedMetrics as PerFlowcellMetrics {
        input:
            aligned_bam    = final_bam,
            aligned_bai    = final_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = metrics_dir
    }

    # Collect stats on decontaminated reads
    if(defined(contaminant_ref_map_file)) {
        # Metrics on contaminated bam
        call SRUTIL.ComputeBamStats as t_020_ComputeContaminatedBamStats { input: bam_file = select_first([DecontaminateSample.contaminated_bam]) }
        # Metrics on original bam file
        call SRUTIL.ComputeBamStats as t_021_ComputeUnalignedBamStats { input: bam_file = select_first([t_002_RevertSam.bam, DecontaminateSample.all_reads_bam])}
    }

    ############################################
    #      _____ _             _ _
    #     |  ___(_)_ __   __ _| (_)_______
    #     | |_  | | '_ \ / _` | | |_  / _ \
    #     |  _| | | | | | (_| | | |/ /  __/
    #     |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ############################################
    File keyfile = t_014_ComputeBamStats.results_file

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
                t_005_AlignReads.bam,
                merged_bam,
                t_007_MarkDuplicates.bam,
                t_008_SortAlignedDuplicateMarkedBam.output_bam,
                t_008_SortAlignedDuplicateMarkedBam.output_bam_index,
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
                t_007_MarkDuplicates.metrics,
                t_011_SamStats.sam_stats,
                t_014_ComputeBamStats.results_file,
                t_015_ComputeBamStatsQ5.results_file,
                t_016_ComputeBamStatsQ7.results_file,
                t_017_ComputeBamStatsQ10.results_file,
                t_018_ComputeBamStatsQ12.results_file,
                t_019_ComputeBamStatsQ15.results_file,
            ],
            keyfile = keyfile
    }

    # Finalize BQSR Metrics if it was run:
    if (perform_BQSR) {
        call FF.FinalizeToDir as t_028_FinalizeBQSRMetrics {
            input:
                outdir = metrics_dir,
                files = select_all([t_009_BaseRecalibrator.recalibration_report]),
                keyfile = keyfile
        }

    }

    call FF.FinalizeToFile as t_029_FinalizeFastQCReport {
        input:
            outdir = metrics_dir,
            file = t_012_FastQC.report
    }

    # Prep a few files for output:
    File fq1_o = unaligned_reads_dir + "/" + basename(fq_e1)
    File fq2_o = unaligned_reads_dir + "/" + basename(fq_e2)
    if (defined(bam)) {
        File unaligned_bam_o = unaligned_reads_dir + "/" + basename(select_first([bam]))
        File unaligned_bai_o = unaligned_reads_dir + "/" + basename(select_first([bai]))
        File fqboup = unaligned_reads_dir + "/" + basename(select_first([DecontaminateSample.decontaminated_unpaired, t_003_Bam2Fastq.fq_unpaired]))
    }
    
    ############################################
    #      ___        _               _
    #     / _ \ _   _| |_ _ __  _   _| |_
    #    | | | | | | | __| '_ \| | | | __|
    #    | |_| | |_| | |_| |_) | |_| | |_
    #     \___/ \__,_|\__| .__/ \__,_|\__|
    #                    |_|
    ############################################
    output {
        # Unaligned reads - If decontamination is ran, these also carry the decontaminated bams
        File fq1 = fq1_o
        File fq2 = fq2_o
        File? fq_unpaired = fqboup

        # Unaligned BAM file
        File? unaligned_bam = unaligned_bam_o
        File? unaligned_bai = unaligned_bai_o

        # Contaminated BAM file and metrics
        File? contaminated_bam = DecontaminateSample.contaminated_bam
        Float? num_contam_reads = select_first([t_020_ComputeContaminatedBamStats.results])['reads']
        Float? pct_contam_reads = select_first([t_020_ComputeContaminatedBamStats.results])['reads'] / select_first([t_021_ComputeUnalignedBamStats.results])['reads'] * 100.0

        # Aligned BAM file
        File aligned_bam = t_025_FinalizeAlignedBam.gcs_path
        File aligned_bai = t_026_FinalizeAlignedBai.gcs_path

        # Unaligned read stats
        Float num_reads = t_014_ComputeBamStats.results['reads']
        Float num_bases = t_014_ComputeBamStats.results['bases']
        Float raw_est_fold_cov = t_014_ComputeBamStats.results['bases']/t_013_ComputeGenomeLength.length

        Float read_length = t_014_ComputeBamStats.results['read_mean']

        Float read_qual_mean = t_014_ComputeBamStats.results['mean_qual']
        Float read_qual_median = t_014_ComputeBamStats.results['median_qual']

        Float num_reads_Q5 = t_015_ComputeBamStatsQ5.results['reads']
        Float num_reads_Q7 = t_016_ComputeBamStatsQ7.results['reads']
        Float num_reads_Q10 = t_017_ComputeBamStatsQ10.results['reads']
        Float num_reads_Q12 = t_018_ComputeBamStatsQ12.results['reads']
        Float num_reads_Q15 = t_019_ComputeBamStatsQ15.results['reads']

        # Aligned read stats
        Float aligned_num_reads = t_012_FastQC.stats_map['number_of_reads']
        Float aligned_num_bases = t_011_SamStats.stats_map['bases_mapped']
        Float aligned_frac_bases = t_011_SamStats.stats_map['bases_mapped']/t_011_SamStats.stats_map['total_length']
        Float aligned_est_fold_cov = t_011_SamStats.stats_map['bases_mapped']/t_013_ComputeGenomeLength.length

        Float aligned_read_length = t_012_FastQC.stats_map['read_length']

        Float insert_size_average = t_011_SamStats.stats_map['insert_size_average']
        Float insert_size_standard_deviation = t_011_SamStats.stats_map['insert_size_standard_deviation']
        Float pct_properly_paired_reads = t_011_SamStats.stats_map['percentage_of_properly_paired_reads_%']

        Float average_identity = 100.0 - (100.0*t_011_SamStats.stats_map['mismatches']/t_011_SamStats.stats_map['bases_mapped'])

        File fastqc_report = t_029_FinalizeFastQCReport.gcs_path
    }
}
