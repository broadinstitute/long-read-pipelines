version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/BAMutils.wdl" as BU

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SplitBamByReadgroup {
    meta {
        description: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
    }
    parameter_meta {
        input_bam: "BAM to operate on; when input is an aligned BAM, user can decide to reset the alignment or not via unmap_bam"
        unmap_bam: "if the input bam should be de-aligned; set to false if input is an uBAM"
        convert_to_fq: "if the bam should be converted to fastq."
        override_workflow_name: "if provided, result files will be saved to gcs_out_root_dir/override_workflow_name; typically you want this when this workflow is called as a subworkflow"
    }
    input {
        File  input_bam
        File? input_bai

        String desired_sample_name

        Boolean unmap_bam
        Boolean convert_to_fq

        Boolean validate_output_bams = false
        Boolean skip_rg_line_pl_value_check = false

        String gcs_out_root_dir
        String? override_workflow_name
        Boolean debug_mode = false
    }
    output {
        Map[String, String] rgid_2_bam = MapRgid2Bams.output_map
        Map[String, Int]    rgid_2_molecule_counts = MapRgid2MoleculeCounts.output_map
        Map[String, String] rgid_2_PU  = MapRgid2PU.output_map
        Map[String, String]? rgid_2_ubam_emptyness = MapRgid2BamEmptiness.output_map
        Boolean rgid_2_bam_are_aligned = ! unmap_bam
        Map[String, String]? rgid_2_fastq = MapRgid2Fqs.output_map

        String last_processing_date = today.yyyy_mm_dd
    }

    ##############################################################################################################################
    String workflow_name = "SplitBamByReadgroup"
    String save_to_this_dir = select_first([override_workflow_name, workflow_name])
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{save_to_this_dir}/" + basename(input_bam, '.bam')
    ##############################################################################################################################

    call AddArtificialReadGroupTag { input: bam = input_bam, desired_sample_name = desired_sample_name }

    ##############################################################################################################################
    # split, if there're > 1 RGs
    call BU.GetReadGroupLines { input: bam = AddArtificialReadGroupTag.bam_with_rg }
    if ( 1 < length(GetReadGroupLines.read_group_ids) ) {
        String output_prefix = basename(AddArtificialReadGroupTag.bam_with_rg, ".bam")
        Int inflation_factor = if(ceil(size([AddArtificialReadGroupTag.bam_with_rg], "GB"))> 150) then 4 else 3
        call Utils.ComputeAllowedLocalSSD as Guess {
            input: intended_gb = 10 + inflation_factor * ceil(size([AddArtificialReadGroupTag.bam_with_rg], "GB"))
        }
        call BU.SplitByRG {
            input:
                bam = AddArtificialReadGroupTag.bam_with_rg, out_prefix = output_prefix, num_ssds = Guess.numb_of_local_ssd
        }
    }
    Array[File] use_these_bams = select_first([SplitByRG.split_bam, [AddArtificialReadGroupTag.bam_with_rg]])

    ##############################################################################################################################
    # unmap/fastq, if so requested
    scatter (bam in use_these_bams) {

        # basecall_model only applies to ONT, so PacBio data will always get 'None'
        Array[String] readgroup_attrs_to_get = ['ID', 'LB', 'PU']
        # call BU.GetReadGroupInfo { input: bam = bam, keys = readgroup_attrs_to_get, null_value_representation = 'None' }
        String rgid = "null" # GetReadGroupInfo.read_group_info['ID']
        String library = "null" # GetReadGroupInfo.read_group_info['LB']
        String platform_unit = "null" # GetReadGroupInfo.read_group_info['PU']

        if (debug_mode) {
            call Utils.CountBamRecords { input: bam = bam }
        }

        # drop alignment if so requested
        if (unmap_bam) {
            call BU.SamtoolsReset as Magic { input: bam = bam }
            call BU.QuerynameSortBamWithPicard as SortUnaligned { input: bam = Magic.res }
            call Utils.CountBamRecords as CountUnalignedRecords { input: bam = Magic.res }
            Boolean uBAM_is_empty = if (0==CountUnalignedRecords.num_records) then true else false

            call FF.FinalizeToFile as SaveUBam {
                input: file = SortUnaligned.qnsort_bam, outdir = outdir
            }
        }
        if (!unmap_bam) {
            call Utils.CountBamRecords as CountInputRGRecords { input: bam = bam }
            Boolean orig_rg_BAM_is_empty = if (0==CountInputRGRecords.num_records) then true else false
            call FF.FinalizeToFile as SaveAlnBam {
                input: file = bam, outdir = outdir
            }
        }
        if (validate_output_bams) {
            Array[String] skip_list_default = ["INVALID_TAG_NM",  # for the purpose we currently have, NM and CIGAR don't matter, and longreads have no mates
                                               "MISSING_TAG_NM",
                                               "INVALID_CIGAR",
                                               "ADJACENT_INDEL_IN_CIGAR",
                                               "CIGAR_MAPS_OFF_REFERENCE",
                                               "MISMATCH_MATE_CIGAR_STRING",
                                               "MATE_CIGAR_STRING_INVALID_PRESENCE",
                                               "MATE_NOT_FOUND",
                                               "INVALID_MAPPING_QUALITY",
                                               "INVALID_FLAG_MATE_UNMAPPED",
                                               "MISMATCH_FLAG_MATE_UNMAPPED",
                                               "INVALID_FLAG_MATE_NEG_STRAND",
                                               "MISMATCH_FLAG_MATE_NEG_STRAND",
                                               "INVALID_MATE_REF_INDEX",
                                               "MISMATCH_MATE_REF_INDEX",
                                               "MISMATCH_MATE_ALIGNMENT_START",
                                               "MATE_FIELD_MISMATCH",
                                               "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND"
                                               ]
            Array[String] for_some_stupid_gcc = ["INVALID_TAG_NM",  # for the purpose we currently have, NM and CIGAR don't matter, and longreads have no mates
                                               "MISSING_TAG_NM",
                                               "INVALID_CIGAR",
                                               "ADJACENT_INDEL_IN_CIGAR",
                                               "CIGAR_MAPS_OFF_REFERENCE",
                                               "MISMATCH_MATE_CIGAR_STRING",
                                               "MATE_CIGAR_STRING_INVALID_PRESENCE",
                                               "MATE_NOT_FOUND",
                                               "INVALID_MAPPING_QUALITY",
                                               "INVALID_FLAG_MATE_UNMAPPED",
                                               "MISMATCH_FLAG_MATE_UNMAPPED",
                                               "INVALID_FLAG_MATE_NEG_STRAND",
                                               "MISMATCH_FLAG_MATE_NEG_STRAND",
                                               "INVALID_MATE_REF_INDEX",
                                               "MISMATCH_MATE_REF_INDEX",
                                               "MISMATCH_MATE_ALIGNMENT_START",
                                               "MATE_FIELD_MISMATCH",
                                               "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND",
                                               "MISSING_PLATFORM_VALUE"
                                               ]
            call BU.ValidateSamFile { input:
                bam = select_first([SortUnaligned.qnsort_bam, bam]),
                validation_errs_to_ignore = if (skip_rg_line_pl_value_check) then skip_list_default else for_some_stupid_gcc
            }
        }

        call BU.CountMolecules as CountMoleculesInRG { input: bam = select_first([SortUnaligned.qnsort_bam, bam]), localize_bam = true }

        # convert to FASTQ if so requested
        if (convert_to_fq) {
            call BU.BamToFastq { input: bam = bam, prefix = basename(bam, ".bam") }
            call FF.FinalizeToFile as SaveFq {
                input: file = BamToFastq.reads_fq, outdir = outdir
            }
            if (debug_mode) { call Utils.CountFastqRecords { input: fastq = BamToFastq.reads_fq } }
        }
    }
    Array[String]  phased_rg_ids   = rgid
    Array[String]  phased_PUs      = platform_unit
    Array[String]  phased_bams     = select_all(flatten([SaveUBam.gcs_path, SaveAlnBam.gcs_path]))
    Array[Boolean] are_ubams_empty = select_all(flatten([uBAM_is_empty, orig_rg_BAM_is_empty]))
    Array[Int]     molecule_counts = CountMoleculesInRG.number

    call GU.CoerceArrayOfPairsToMap as MapRgid2PU { input: keys = phased_rg_ids, values = phased_PUs }
    call GU.CoerceArrayOfPairsToMap as MapRgid2Bams { input: keys = phased_rg_ids, values = phased_bams }
    if (convert_to_fq) {
        Array[String] phased_fastqs = select_all(SaveFq.gcs_path)
        call GU.CoerceArrayOfPairsToMap as MapRgid2Fqs { input: keys = phased_rg_ids, values = phased_fastqs }
    }
    if (unmap_bam) {
        call GU.CoerceArrayOfPairsToMap as MapRgid2BamEmptiness { input: keys = phased_rg_ids, values = are_ubams_empty }
    }
    call GU.CoerceArrayOfPairsToMap as MapRgid2MoleculeCounts { input: keys = phased_rg_ids, values = molecule_counts }

    call GU.GetTodayDate as today {}
}

task AddArtificialReadGroupTag {
    input {
        File bam
        String desired_sample_name
    }
    output {
        File bam_with_rg = "~{prefix}_withArtificialRG.bam"
        File bai_with_rg = "~{prefix}_withArtificialRG.bam.bai"
    }
    String prefix = basename(bam, ".bam")
    command <<<
    set -euxo pipefail

        # # handle the header
        # echo -e "@RG\tID:artificial\tSM:~{desired_sample_name}\tPL:ONT\tPU:null\tLB:null" \
        # > use_this_rg_line.txt
        # samtools view -H ~{bam} \
        # > header.txt
        # cat header.txt \
        #     use_this_rg_line.txt \
        # > new_header.txt

        # add RG tag to each read, and reheader at the same time to add the RG header line
        samtools view -@3 -h "~{bam}" \
        | awk 'BEGIN {
            OFS="\t"
            rg="@RG\tID:artificial\tSM:~{desired_sample_name}\tPL:ONT\tPU:null\tLB:null"
            }
            /^@/ { print; next }
            !rg_printed { print rg; rg_printed=1 }
            { print $0, "RG:Z:artificial" }' \
        | samtools view -@3 --write-index \
            -u \
            -o "~{prefix}_withArtificialRG.bam##idx##~{prefix}_withArtificialRG.bam.bai" \
            -

        # verify
        samtools view -H "~{prefix}_withArtificialRG.bam" \
        | grep -c "^@RG"
        samtools view -H "~{prefix}_withArtificialRG.bam" \
        | grep "^@RG"

        samtools view \
            "~{prefix}_withArtificialRG.bam" \
        | head \
        | grep -cF "RG:Z:artificial" || true
    >>>
    runtime {
        cpu:            6
        memory:         "24 GiB"
        disks:          "local-disk 750 LOCAL"
        preemptible:    0
        maxRetries:     0
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23.1"
    }
}
