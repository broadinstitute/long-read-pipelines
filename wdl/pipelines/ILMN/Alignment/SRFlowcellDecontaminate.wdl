version 1.0

import "../../../tasks/Utility/SRUtils.wdl" as SRUTIL
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Preprocessing/RemoveSingleOrganismContamination.wdl" as DECONTAMINATE
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRFlowcellDecontaminate {

    meta {
        author: "Raphael Brosula, script modified form Jonn Smith"
        description: "This workflow preprocesses BAMs or paired fastq files for decontamination ONLY."
    }

    parameter_meta {
        bam: "Bam file containing reads to align and process.  `bai` must be defined if this argument is.  This argument and `bai` are mutually exclusive with `fq_end1` and `fq_end2`"
        bai: "Index for `bam`.  `bam` must be defined if this argument is.  This argument and `bam` are mutually exclusive with `fq_end1` and `fq_end2`"

        fq_end1: "FASTQ file containing end 1 of the short read data to process.  `fq_end2` must be defined if this argument is.  This argument and `fq_end2` are mutually exclusive with `bam` and `bai`"
        fq_end2: "FASTQ file containing end 2 of the short read data to process.  `fq_end1` must be defined if this argument is.  This argument and `fq_end1` are mutually exclusive with `bam` and `bai`"

        SM: "Sample name for the given bam file."
        LB: "Library name for the given bam file."

        contaminant_ref_name:  "Name for the contaminant reference."
        contaminant_ref_map_file:  "Reference map file for the contaminant reference sequence and auxillary file locations."

        dir_prefix: "Directory prefix to use for finalized location."
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs."

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

        String? contaminant_ref_name
        File? contaminant_ref_map_file

        String dir_prefix

        String gcs_out_root_dir

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

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    # Create an outdir:
    String outdir = if DEBUG_MODE then sub(gcs_out_root_dir, "/$", "") + "/SRFlowcell/~{dir_prefix}/" + t_001_WdlExecutionStartTimestamp.timestamp_string else sub(gcs_out_root_dir, "/$", "") + "/SRFlowcell/~{dir_prefix}"
    String reads_dir = outdir + "/reads"

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
        call DECONTAMINATE.RemoveSingleOrganismContamination as t_005_DecontaminateSample {
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

    File fq1 = select_first([t_005_DecontaminateSample.decontaminated_fq1, fq_end1, t_003_Bam2Fastq.fq_end1])
    File fq2 = select_first([t_005_DecontaminateSample.decontaminated_fq2, fq_end2, t_003_Bam2Fastq.fq_end2])
    File contaminated_bam = select_first([t_005_DecontaminateSample.contaminated_bam])

    # Metrics on contaminated_bam
    call SRUTIL.ComputeBamStats as t_006_ComputeContaminatedBamStats { input: bam_file = contaminated_bam }
    # Metrics on original files
    call SRUTIL.ComputeBamStats as t_007_ComputeUnalignedBamStats { input: bam_file = select_first([t_002_RevertSam.bam, t_005_DecontaminateSample.all_reads_bam])}

    ############################################
    #      _____ _             _ _
    #     |  ___(_)_ __   __ _| (_)_______
    #     | |_  | | '_ \ / _` | | |_  / _ \
    #     |  _| | | | | | (_| | | |/ /  __/
    #     |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ############################################
    call FF.FinalizeToFile as t_008_FinalizeFq1 { input: outdir = reads_dir, file = fq1, keyfile = contaminated_bam }
    call FF.FinalizeToFile as t_009_FinalizeFq2 { input: outdir = reads_dir, file = fq2, keyfile = contaminated_bam }
    call FF.FinalizeToFile as t_010_FinalizeContaminatedBam{ input: outdir = reads_dir, file = contaminated_bam, keyfile = contaminated_bam }

    output {
        File fq_o1 = t_008_FinalizeFq1.gcs_path
        File fq_o2 = t_009_FinalizeFq2.gcs_path
        File contaminated_bamout = t_010_FinalizeContaminatedBam.gcs_path

        Float num_total_reads = t_007_ComputeUnalignedBamStats.results['reads']
        Float num_contam_reads = t_006_ComputeContaminatedBamStats.results['reads']
        Float pct_contam_reads = t_006_ComputeContaminatedBamStats.results['reads'] / t_007_ComputeUnalignedBamStats.results['reads'] * 100.0
    }
}
