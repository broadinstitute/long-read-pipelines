version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/BAMutils.wdl" as BU

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SplitBamByReadgroup {
    meta {
        desciption: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
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

        Boolean unmap_bam
        Boolean convert_to_fq

        Boolean validate_output_bams = false

        String gcs_out_root_dir
        String? override_workflow_name
        Boolean debug_mode = false
    }
    output {
        Map[String, String] rgid_2_bam = MapRgid2Bams.output_map
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

    ##############################################################################################################################
    # split, if there're > 1 RGs
    call BU.GetReadGroupLines { input: bam = input_bam }
    if ( 1 < length(GetReadGroupLines.read_group_ids) ) {
        String output_prefix = basename(input_bam, ".bam")
        Int inflation_factor = if(ceil(size([input_bam], "GB"))> 150) then 4 else 3
        call Utils.ComputeAllowedLocalSSD as Guess {
            input: intended_gb = 10 + inflation_factor * ceil(size([input_bam], "GB"))
        }
        call BU.SplitByRG {
            input:
                bam = input_bam, out_prefix = output_prefix, num_ssds = Guess.numb_of_local_ssd
        }
    }
    Array[File] use_these_bams = select_first([SplitByRG.split_bam, [input_bam]])

    ##############################################################################################################################
    # unmap/fastq, if so requested
    scatter (bam in use_these_bams) {

        # basecall_model only applies to ONT, so PacBio data will always get 'None'
        Array[String] readgroup_attrs_to_get = ['ID', 'LB', 'PU']
        call BU.GetReadGroupInfo { input: bam = bam, keys = readgroup_attrs_to_get, null_value_representation = 'None' }
        String rgid = GetReadGroupInfo.read_group_info['ID']
        String library = GetReadGroupInfo.read_group_info['LB']
        String platform_unit = GetReadGroupInfo.read_group_info['PU']

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
            call BU.ValidateSamFile { input: bam = select_first([SortUnaligned.qnsort_bam, bam]) }
        }

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

    call GU.CoerceArrayOfPairsToMap as MapRgid2PU { input: keys = phased_rg_ids, values = phased_PUs }
    call GU.CoerceArrayOfPairsToMap as MapRgid2Bams { input: keys = phased_rg_ids, values = phased_bams }
    if (convert_to_fq) {
        Array[String] phased_fastqs = select_all(SaveFq.gcs_path)
        call GU.CoerceArrayOfPairsToMap as MapRgid2Fqs { input: keys = phased_rg_ids, values = phased_fastqs }
    }
    if (unmap_bam) {
        call GU.CoerceArrayOfPairsToMap as MapRgid2BamEmptiness { input: keys = phased_rg_ids, values = are_ubams_empty }
    }

    call GU.GetTodayDate as today {}
}
