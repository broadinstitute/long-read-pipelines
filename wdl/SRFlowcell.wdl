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
import "tasks/NanoPlot.wdl" as NP
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

        String dir_prefix

        String gcs_out_root_dir

        Boolean DEBUG_MODE = false
    }

    parameter_meta {
        bam:                "GCS path to unmapped bam"
        bai:                "GCS path to bai index for unmapped bam"

        fq_end1:            "GCS path to end1 of paired-end fastq"
        fq_end2:            "GCS path to end2 of paired-end fastq"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"

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
    }

    File fq1 = select_first([fq_end1, t_003_Bam2Fastq.fq_end1])
    File fq2 = select_first([fq_end2, t_003_Bam2Fastq.fq_end2])

    # Align reads to reference with BWA-MEM2:
    call SRUTIL.BwaMem2 as t_004_AlignReads {
        input:
            fq_end1 = fq1,
            fq_end2 = fq2,
            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],
            ref_0123 = ref_map["0123"],
            ref_amb = ref_map["amb"],
            ref_ann = ref_map["ann"],
            ref_bwt = ref_map["bwt"],
            ref_pac = ref_map["pac"],
            mark_short_splits_as_secondary = true,
            prefix = SM + ".aligned"
    }

    if (defined(bam)) {
        # Merge aligned reads and unaligned reads:
        call SRUTIL.MergeBamAlignment as t_005_MergeBamAlignment {
            input:
                aligned_bam = t_004_AlignReads.bam,
                unaligned_bam = select_first([t_002_RevertSam.bam]),
                ref_fasta = ref_map["fasta"],
                ref_fasta_index = ref_map["fai"],
                ref_dict = ref_map["dict"],
                prefix = SM + ".aligned.merged"
        }
    }

    File merged_bam = select_first([t_004_AlignReads.bam, t_005_MergeBamAlignment.bam])

# This was for GATK 3.  We probably don't need it now.
#    # Fix mates:
#    call SRUTIL.FixMate as t_006_FixMate {
#        input:
#            input_bam = merged_bam,
#            prefix = SM + ".aligned.merged.fixmates"
#    }

    # Mark Duplicates
    call SRUTIL.MarkDuplicates as t_006_MarkDuplicates {
        input:
            input_bam = merged_bam,
            prefix = SM + ".aligned.merged.markDuplicates"
    }

    # Sort Duplicate Marked Bam:
    call Utils.SortBam as t_007_SortAlignedDuplicateMarkedBam {
        input:
            input_bam = t_006_MarkDuplicates.bam,
            prefix = SM + ".aligned.merged.markDuplicates.sorted"
    }

#    Fingerprinting?

    # Recalibrate Base Scores:
    call SRUTIL.BaseRecalibrator as t_008_BaseRecalibrator {
        input:
            input_bam = t_007_SortAlignedDuplicateMarkedBam.sorted_bam,
            input_bam_index = t_007_SortAlignedDuplicateMarkedBam.sorted_bai,

            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],

            known_sites_vcf = ref_map["known_sites_vcf"],
            known_sites_index = ref_map["known_sites_index"],

            prefix = SM + ".baseRecalibratorReport"
    }

    call SRUTIL.ApplyBQSR as t_009_ApplyBQSR {
        input:
            input_bam = t_007_SortAlignedDuplicateMarkedBam.sorted_bam,
            input_bam_index = t_007_SortAlignedDuplicateMarkedBam.sorted_bai,

            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],

            recalibration_report = t_008_BaseRecalibrator.recalibration_report,

            prefix = SM + ".aligned.merged.markDuplicates.sorted.BQSR"
    }

    #############################################
    #      __  __      _        _
    #     |  \/  | ___| |_ _ __(_) ___ ___
    #     | |\/| |/ _ \ __| '__| |/ __/ __|
    #     | |  | |  __/ |_| |  | | (__\__ \
    #     |_|  |_|\___|\__|_|  |_|\___|___/
    #
    #############################################

    call AM.SamStats as t_010_SamStats {
        input:
            bam = t_009_ApplyBQSR.recalibrated_bam
    }

    call NP.NanoPlotFromBam as t_011_NanoPlotFromBam {
        input:
            bam = t_009_ApplyBQSR.recalibrated_bam,
            bai = t_009_ApplyBQSR.recalibrated_bai
    }

    call Utils.ComputeGenomeLength as t_012_ComputeGenomeLength {
        input:
            fasta = ref_map['fasta']
    }

    ############################################
    #      _____ _             _ _
    #     |  ___(_)_ __   __ _| (_)_______
    #     | |_  | | '_ \ / _` | | |_  / _ \
    #     |  _| | | | | | (_| | | |/ /  __/
    #     |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ############################################
    File keyfile = t_010_SamStats.sam_stats
    String reads_dir = outdir + "/reads"
    String metrics_dir = outdir + "/metrics"

    # Finalize our reads first:
#    call FF.FinalizeToDir as t_013_FinalizeUnalignedReads {
#        input:
#            outdir = reads_dir + "/unaligned",
#            files =
#            [
#                bam,
#                bai,
#                t_003_Bam2Fastq.fastq
#            ],
#            keyfile = keyfile
#    }

    call FF.FinalizeToDir as t_014_FinalizeAlignedReads {
        input:
            outdir = reads_dir + "/aligned",
            files =
            [
                t_004_AlignReads.bam,
                merged_bam,
                t_006_MarkDuplicates.bam,
                t_007_SortAlignedDuplicateMarkedBam.sorted_bam,
                t_007_SortAlignedDuplicateMarkedBam.sorted_bai,
            ],
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_015_FinalizeAlignedBam {
        input:
            outdir = reads_dir + "/aligned",
            file = t_009_ApplyBQSR.recalibrated_bam,
            keyfile = keyfile
    }

    call FF.FinalizeToFile as t_016_FinalizeAlignedBai {
        input:
            outdir = reads_dir + "/aligned",
            file = t_009_ApplyBQSR.recalibrated_bai,
            keyfile = keyfile
    }

    # Finalize our metrics:
    call FF.FinalizeToDir as t_017_FinalizeMetrics {
        input:
            outdir = metrics_dir,
            files =
            [
                t_006_MarkDuplicates.metrics,
                t_008_BaseRecalibrator.recalibration_report,
                t_010_SamStats.sam_stats
            ],
            keyfile = keyfile
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
        # Aligned BAM file
        File aligned_bam = t_015_FinalizeAlignedBam.gcs_path
        File aligned_bai = t_016_FinalizeAlignedBai.gcs_path

        # Aligned read stats
        Float aligned_num_reads = t_011_NanoPlotFromBam.stats_map['number_of_reads']
        Float aligned_num_bases = t_011_NanoPlotFromBam.stats_map['number_of_bases_aligned']
        Float aligned_frac_bases = t_011_NanoPlotFromBam.stats_map['fraction_bases_aligned']
        Float aligned_est_fold_cov = t_011_NanoPlotFromBam.stats_map['number_of_bases_aligned']/t_012_ComputeGenomeLength.length

        Float aligned_read_length_mean = t_011_NanoPlotFromBam.stats_map['mean_read_length']
        Float aligned_read_length_median = t_011_NanoPlotFromBam.stats_map['median_read_length']
        Float aligned_read_length_stdev = t_011_NanoPlotFromBam.stats_map['read_length_stdev']
        Float aligned_read_length_N50 = t_011_NanoPlotFromBam.stats_map['n50']

        Float average_identity = t_011_NanoPlotFromBam.stats_map['average_identity']
        Float median_identity = t_011_NanoPlotFromBam.stats_map['median_identity']
    }
}
