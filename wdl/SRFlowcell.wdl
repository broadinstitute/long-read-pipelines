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
import "tasks/Finalize.wdl" as FF

workflow SRFlowcell {
    input {
        File bam
        File bai

        String SM
        String LB

        File ref_map_file

        Int? num_shards
        String dir_prefix

        String gcs_out_root_dir

        Boolean DEBUG_MODE = false
    }

    parameter_meta {
        bam:                "GCS path to unmapped bam"
        bai:                "GCS path to bai index for unmapped bam"
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
    call Utils.GetCurrentTimestampString as WdlExecutionStartTimestamp { input: }

    # Create an outdir:
    String outdir = if DEBUG_MODE then sub(gcs_out_root_dir, "/$", "") + "/SRFlowcell/~{dir_prefix}/" + WdlExecutionStartTimestamp.timestamp_string else sub(gcs_out_root_dir, "/$", "") + "/SRFlowcell/~{dir_prefix}"

    # Convert input SAM/BAM to FASTQ:
    call SRUTIL.Bam2FqPicard as Bam2Fastq {
        input:
            bam = bam,
            prefix = SM
    }

    # Align reads to reference with BWA-MEM2:
    call SRUTIL.BwaMem2 as AlignReads {
        input:
            name_sorted_fastq = Bam2Fastq.fastq,
            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],
            ref_0123 = ref_map["0123"],
            ref_amb = ref_map["amb"],
            ref_ann = ref_map["ann"],
            ref_bwt = ref_map["bwt"],
            ref_pac = ref_map["pac"],
            prefix = SM + ".aligned"
    }

    # Merge aligned reads and unaligned reads:
    call SRUTIL.MergeBamAlignment as MergeBamAlignment {
        input:
            aligned_bam = AlignReads.bam,
            unaligned_bam = bam,
            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],
    }

    # Mark Duplicates
    call SRUTIL.MarkDuplicates as MarkDuplicates {
        input:
            input_bam = MergeBamAlignment.bam,
            prefix = SM + ".aligned.markDuplicates."
    }

    # Sort Duplicate Marked Bam:
    call Utils.SortBam as SortAlignedDuplicateMarkedBam {
        input:
            input_bam = MarkDuplicates.bam,
            prefix = SM + ".aligned.markDuplicates.sorted"
    }

#    Fingerprinting?
#
    # BQSR
    call SRUTIL.BaseRecalibrator as BaseRecalibrator {
        input:
            input_bam = SortAlignedDuplicateMarkedBam.sorted_bam,
            input_bam_index = SortAlignedDuplicateMarkedBam.sorted_bai,

            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],

            known_sites_vcf = ref_map["known_sites_vcf"],
            known_sites_index = ref_map["known_sites_vcf_idx"],

            prefix = SM + ".baseRecalibratorReport"
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

    }
}