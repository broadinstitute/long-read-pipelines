version 1.0

##########################################################################################
## A workflow that preprocesses short read flowcell data in preparation for variant calling.
## This workflow contains the following steps:
##   1) Sam -> Fastq (if necessary)
##   2) Alignment to reference with bwa-mem2 (https://github.com/bwa-mem2/bwa-mem2)
##   3) Mark Duplicate reads
##   4) Recalibrate base quality scores.
##########################################################################################

import "tasks/SRUtils.wdl" as SRUTIL
import "tasks/Utils.wdl" as Utils
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/FastQC.wdl" as FastQC
import "tasks/RemoveSingleOrganismContamination.wdl" as DECONTAMINATE
import "tasks/Finalize.wdl" as FF

workflow SRFlowcell {
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

    parameter_meta {
        bam:                "GCS path to unmapped bam"
        bai:                "GCS path to bai index for unmapped bam"

        fq_end1:            "GCS path to end1 of paired-end fastq"
        fq_end2:            "GCS path to end2 of paired-end fastq"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        contaminant_ref_name:                 "Name of the contaminant genome to be used in output files."
        contaminant_ref_map_file:       "table indicating reference sequence and auxillary file locations for a single-organism contaminant"

        SM:                 "the value to place in the BAM read group's SM field"
        LB:                 "the value to place in the BAM read group's LB (library) field"

        num_shards:         "number of shards into which fastq files should be batched"
        dir_prefix:         "directory prefix for output files"

        DEBUG_MODE:         "[default valued] enables debugging tasks / subworkflows (default: false)"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
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
    if (false && defined(contaminant_ref_map_file)) {

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
                input_bam = t_008_SortAlignedDuplicateMarkedBam.sorted_bam,
                input_bam_index = t_008_SortAlignedDuplicateMarkedBam.sorted_bai,

                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],

                known_sites_vcf = ref_map["known_sites_vcf"],
                known_sites_index = ref_map["known_sites_index"],

                prefix = SM + ".baseRecalibratorReport"
        }

        call SRUTIL.ApplyBQSR as t_010_ApplyBQSR {
            input:
                input_bam = t_008_SortAlignedDuplicateMarkedBam.sorted_bam,
                input_bam_index = t_008_SortAlignedDuplicateMarkedBam.sorted_bai,

                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],

                recalibration_report = t_009_BaseRecalibrator.recalibration_report,

                prefix = SM + ".aligned.merged.markDuplicates.sorted.BQSR"
        }
    }

    File final_bam = select_first([t_010_ApplyBQSR.recalibrated_bam, t_008_SortAlignedDuplicateMarkedBam.sorted_bam])
    File final_bai = select_first([t_010_ApplyBQSR.recalibrated_bai, t_008_SortAlignedDuplicateMarkedBam.sorted_bai])

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
    call FF.FinalizeToDir as t_020_FinalizeUnalignedFastqReads {
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
        call FF.FinalizeToDir as t_021_FinalizeUnalignedReadsFromBam {
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

    call FF.FinalizeToDir as t_022_FinalizeAlignedReads {
        input:
            outdir = aligned_reads_dir,
            files =
            [
                t_005_AlignReads.bam,
                merged_bam,
                t_007_MarkDuplicates.bam,
                t_008_SortAlignedDuplicateMarkedBam.sorted_bam,
                t_008_SortAlignedDuplicateMarkedBam.sorted_bai,
            ],
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_023_FinalizeAlignedBam {
        input:
            outdir = aligned_reads_dir,
            file = final_bam,
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_024_FinalizeAlignedBai {
        input:
            outdir = aligned_reads_dir,
            file = final_bai,
            keyfile = keyfile
    }

    # Finalize our metrics:
    call FF.FinalizeToDir as t_025_FinalizeMetrics {
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
        call FF.FinalizeToDir as t_026_FinalizeBQSRMetrics {
            input:
                outdir = metrics_dir,
                files = select_all([t_009_BaseRecalibrator.recalibration_report]),
                keyfile = keyfile
        }

    }

    call FF.FinalizeToFile as t_027_FinalizeFastQCReport {
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
        # Unaligned reads
        File fq1 = fq1_o
        File fq2 = fq2_o
        File? fq_unpaired = fqboup

        # Unaligned BAM file
        File? unaligned_bam = unaligned_bam_o
        File? unaligned_bai = unaligned_bai_o

        # Contaminated BAM file:
        File? contaminated_bam = DecontaminateSample.contaminated_bam

        # Aligned BAM file
        File aligned_bam = t_023_FinalizeAlignedBam.gcs_path
        File aligned_bai = t_024_FinalizeAlignedBai.gcs_path

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

        File fastqc_report = t_027_FinalizeFastQCReport.gcs_path
    }
}
