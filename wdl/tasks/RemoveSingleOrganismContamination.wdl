version 1.0

import "Structs.wdl" as Structs
import "SRUtils.wdl" as SRUTIL
import "Utils.wdl" as Utils
import "Finalize.wdl" as FF

workflow RemoveSingleOrganismContamination {
    meta {
        author: "Jonn Smith"
        description: "A workflow to remove contamination originating from a single organism from a dataset."
    }

    input {
        File? input_bam
        File? input_bai

        File? fq_end1
        File? fq_end2

        String SM
        String LB
        String platform = "illumina"

        File contaminant_ref_name
        File contaminant_ref_map_file

        String dir_prefix
        String gcs_out_root_dir

        Boolean DEBUG_MODE = false
    }

    parameter_meta {
        input_bam:                "GCS path to unmapped bam"
        input_bai:                "GCS path to bai index for unmapped bam"

        fq_end1:            "GCS path to end1 of paired-end fastq"
        fq_end2:            "GCS path to end2 of paired-end fastq"

        SM:                 "the value to place in the BAM read group's SM field"
        LB:                 "the value to place in the BAM read group's LB (library) field"
        platform:                 "[default valued] the value to place in the BAM read group's PL (platform) field (default: illumina)"

        contaminant_ref_name:                 "Name of the contaminant genome to be used in output files."
        contaminant_ref_map_file:                 "Table indicating reference sequence and auxillary file locations."

        dir_prefix:         "directory prefix for output files"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        DEBUG_MODE:         "[default valued] enables debugging tasks / subworkflows (default: false)"
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    # Some basic error handling:
    if (!defined(input_bam) && (!defined(fq_end1) || !defined(fq_end2))) {
        call Utils.FailWithWarning as t_002_NoInputFileProvidedFailure {
            input: warning = "No input file has been provided!  You must provide either an input bam or input fastq1/fastq2 files."
        }
    }
    if (defined(input_bam) && (defined(fq_end1) || defined(fq_end2))) {
        call Utils.FailWithWarning as t_002_TooManyInputsProvidedFailure {
            input: warning = "Too many inputs provided!  You must provide EITHER an input bam OR input fastq1/fastq2 files."
        }
    }

    # Create an outdir:
    String outdir = if DEBUG_MODE then sub(gcs_out_root_dir, "/$", "") + "/RemoveSingleOrganismContamination/~{dir_prefix}/" + t_001_WdlExecutionStartTimestamp.timestamp_string else sub(gcs_out_root_dir, "/$", "") + "/RemoveSingleOrganismContamination/~{dir_prefix}"

    # Get ref info:
    Map[String, String] ref_map = read_map(contaminant_ref_map_file)

    if (defined(input_bam)) {
        # Convert the given bam to a uBAM (needed for previous aligned data):
        call SRUTIL.RevertSam as t_003_RevertSam {
            input:
                input_bam = select_first([input_bam]),
                prefix = SM + ".revertSam"
        }

        # Convert input SAM/BAM to FASTQ:
        call SRUTIL.BamToFq as t_004_Bam2Fastq {
            input:
                bam = t_003_RevertSam.bam,
                prefix = SM
        }
        call Utils.GetRawReadGroup as t_005_GetRawReadGroup { input: gcs_bam_path = select_first([input_bam]) }
    }

    File fq_e1 = select_first([fq_end1, t_004_Bam2Fastq.fq_end1])
    File fq_e2 = select_first([fq_end2, t_004_Bam2Fastq.fq_end2])

    String RG = select_first([t_005_GetRawReadGroup.rg, "@RG\tID:" + SM + "_" + LB + "\tPL:" + platform + "\tLB:" + LB + "\tSM:" + SM])

    # Align data to contaminant reference:
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
            prefix = SM + ".contaminant_aligned." + contaminant_ref_name
    }

    call Utils.FilterReadsBySamFlags as t_007_ExtractDecontaminatedReads {
        input:
            bam = t_006_AlignReads.bam,
            sam_flags = "256",
            extra_args = " -f 12 ",
            prefix = SM + ".decontaminated"
    }

    call Utils.FilterReadsBySamFlags as t_008_ExtractContaminatedReads {
        input:
            bam = t_006_AlignReads.bam,
            sam_flags = "12",
            prefix = SM + ".contaminated_" + contaminant_ref_name + "_reads"
    }

    call Utils.SortBam as t_009_SortDecontaminatedReads {
        input:
            input_bam = t_007_ExtractDecontaminatedReads.output_bam,
            extra_args = " -n ",
            prefix = SM + ".decontaminated.sorted"
    }

    call Utils.SortBam as t_010_SortContaminatedReads {
        input:
            input_bam = t_008_ExtractContaminatedReads.output_bam,
            extra_args = " -n ",
            prefix = SM + ".contaminated_" + contaminant_ref_name + "_reads.sorted"
    }

    # Convert input SAM/BAM to FASTQ:
    call SRUTIL.BamToFq as t_011_CreateFastqFromDecontaminatedReads {
        input:
            bam = t_009_SortDecontaminatedReads.sorted_bam,
            prefix = SM + ".decontaminated"
    }

    ############################################
    #      _____ _             _ _
    #     |  ___(_)_ __   __ _| (_)_______
    #     | |_  | | '_ \ / _` | | |_  / _ \
    #     |  _| | | | | | (_| | | |/ /  __/
    #     |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ############################################

    # Chosen because it's a relatively small file.
    File keyfile = t_011_CreateFastqFromDecontaminatedReads.monitoring_log

    call FF.FinalizeToFile as t_012_FinalizeContaminatedBam { input: outdir = outdir, file = t_010_SortContaminatedReads.sorted_bam, keyfile = keyfile }
    call FF.FinalizeToFile as t_013_FinalizeContaminatedBai { input: outdir = outdir, file = t_010_SortContaminatedReads.sorted_bai, keyfile = keyfile }
    call FF.FinalizeToFile as t_014_FinalizeDecontaminatedFq1 { input: outdir = outdir, file = t_011_CreateFastqFromDecontaminatedReads.fq_end1, keyfile = keyfile }
    call FF.FinalizeToFile as t_015_FinalizeDecontaminatedFq2 { input: outdir = outdir, file = t_011_CreateFastqFromDecontaminatedReads.fq_end2, keyfile = keyfile }

    ############################################
    #      ___        _               _
    #     / _ \ _   _| |_ _ __  _   _| |_
    #    | | | | | | | __| '_ \| | | | __|
    #    | |_| | |_| | |_| |_) | |_| | |_
    #     \___/ \__,_|\__| .__/ \__,_|\__|
    #                    |_|
    ############################################

    output {
        File contaminated_bam = t_012_FinalizeContaminatedBam.gcs_path
        File contaminated_bam_index = t_013_FinalizeContaminatedBai.gcs_path

        File decontaminated_fq1 = t_014_FinalizeDecontaminatedFq1.gcs_path
        File decontaminated_fq2 = t_015_FinalizeDecontaminatedFq2.gcs_path
    }
}