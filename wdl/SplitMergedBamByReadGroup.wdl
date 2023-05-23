version 1.0

import "tasks/Utils.wdl"
import "tasks/utils/BAMutils.wdl" as BU
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/utils/UnAlignBam.wdl" as UnAlign
import "tasks/Finalize.wdl" as FF
import "tasks/ONTUtils.wdl" as OU

workflow SplitMergedBamByReadGroup {
    meta {
        description: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
    }

    input {
        File input_bam
        File? input_bai

        String platform
        Boolean unmap_bam
        Boolean convert_to_fq = false
        Boolean fix_bam_header = false

        String gcs_out_root_dir

        Boolean debug_mode = false
    }

    parameter_meta {
        input_bam: "BAM to be split by read group; doesn't necessarily need to be aligned or sorted."
        input_bai: "(optional) BAI accompanying the BAM"
        platform:  "long reads platform the BAM was generated on; must be one of [PB, ONT]"
        convert_to_fq: "user option to convert to FASTQ (gz) or not"
        gcs_out_root_dir: "place to store the result files"
    }

    if (platform!="PB" && platform!="ONT") {
        String formatted_reason = "Provided value for 'platform' (" + platform + ") isn't supported. Must be one of [PB, ONT]."
        call Utils.StopWorkflow { input: reason = formatted_reason }
    }

    String workflow_name = "SplitMergedBamByReadGroup"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/" + basename(input_bam, '.bam')

    if (fix_bam_header) {
        call FixParticularBamHeaderIssue { input: bam = input_bam }
    }
    File ok_input_bam = select_first([FixParticularBamHeaderIssue.fixed, input_bam])

    # most basic metadata and QC check
    call BU.GatherBamMetadata {
        input: bam = ok_input_bam
    }
    call Utils.InferSampleName { input: bam = ok_input_bam, bai = input_bai }
    call BU.ValidateSamFile { input: bam = ok_input_bam }

    if ('PB' == platform) {
        # this guarantees that there are no records without RG tag
        call BU.VerifyPacBioBamHasAppropriatePrimroseRuns as PrimroseCheck { input: bam = ok_input_bam }
        if (0!= length(PrimroseCheck.readgroups_missing_primrose)) {
            call Utils.StopWorkflow as MissingPrimrose { input: reason = "Input BAM file has some of its read groups missing primrose calls."}
        }
    }

    # split, if there's > 1 RGs
    call BU.GetReadGroupLines { input: bam = ok_input_bam }
    if ( 1 < length(GetReadGroupLines.read_group_ids) ) {
        String output_prefix = basename(ok_input_bam, ".bam")
        Int inflation_factor = if(ceil(size([ok_input_bam], "GB"))> 150) then 4 else 3
        call Utils.ComputeAllowedLocalSSD as Guess {
            input: intended_gb = 10 + inflation_factor * ceil(size([ok_input_bam], "GB"))
        }
        call BU.SplitByRG {
            input:
                bam = ok_input_bam, out_prefix = output_prefix, num_ssds = Guess.numb_of_local_ssd
        }
    }
    Array[File] use_these_bams = select_first([SplitByRG.split_bam, [ok_input_bam]])

    scatter (bam in use_these_bams) {

        # basecall_model only applies to ONT, so PacBio data will always get 'None'
        Array[String] readgroup_attrs_to_get = ['ID', 'LB', 'PU']
        call BU.GetReadGroupInfo { input: uBAM = bam, keys = readgroup_attrs_to_get, null_value_representation = 'None' }
        String rgid = GetReadGroupInfo.read_group_info['ID']
        String library = GetReadGroupInfo.read_group_info['LB']
        String platform_unit = GetReadGroupInfo.read_group_info['PU']

        if (debug_mode) {
            call Utils.CountBamRecords { input: bam = bam }
        }

        # drop alignment if so requested
        if (unmap_bam) {
            call BU.SamtoolsReset as Magic { input: bam = bam }
            call BU.SortBamByQName as SortUnaligned { input: bam = Magic.res }
            call Utils.CountBamRecords as CountUnalignedRecords { input: bam = Magic.res }
            Boolean uBAM_is_empty = if (0==CountUnalignedRecords.num_records) then true else false

            call FF.FinalizeToFile as SaveUBam {
                input: file = SortUnaligned.qnsort_bam, outdir = outdir
            }
        }
        if (!unmap_bam) {
            call FF.FinalizeToFile as SaveAlnBam {
                input: file = bam, outdir = outdir
            }
        }

        # convert to FASTQ if so requested
        if (convert_to_fq) {
            call Utils.BamToFastq { input: bam = bam, prefix = basename(bam, ".bam") }
            call FF.FinalizeToFile as SaveFq {
                input: file = BamToFastq.reads_fq, outdir = outdir
            }
            if (debug_mode) { call Utils.CountFastqRecords { input: fastq = BamToFastq.reads_fq } }
        }
    }
    Array[String]  phased_rg_ids   = rgid
    Array[String]  phased_PUs      = platform_unit
    Array[String]  phased_bams     = select_first([select_all(SaveUBam.gcs_path), select_all(SaveAlnBam.gcs_path)])
    Array[Boolean?] are_ubams_empty = uBAM_is_empty
    Array[String]? phased_fastqs   = select_first([select_all(SaveFq.gcs_path), select_all(SaveFq.gcs_path)])

    call GU.CoerceArrayOfPairsToMap as MapRgid2PU { input: keys = phased_rg_ids, values = phased_PUs }
    call GU.CoerceArrayOfPairsToMap as MapRgid2Bams { input: keys = phased_rg_ids, values = phased_bams }
    if (convert_to_fq) {
        call GU.CoerceArrayOfPairsToMap as MapRgid2Fqs { input: keys = phased_rg_ids, values = select_first([phased_fastqs]) }
    }
    if (unmap_bam) {
        call GU.CoerceArrayOfPairsToMap as MapRgid2BamEmptiness { input: keys = phased_rg_ids, values = select_all(are_ubams_empty) }
    }

    if (platform=="ONT") {
        call OU.GetBasecallModel { input: bam = ok_input_bam }
    }

    call GU.GetTodayDate as today {}

    output {
        Map[String, String] rgid_2_bam = MapRgid2Bams.output_map
        Map[String, String] rgid_2_PU  = MapRgid2PU.output_map
        Map[String, String]? rgid_2_ubam_emptyness = MapRgid2BamEmptiness.output_map
        Boolean rgid_2_bam_are_aligned = ! unmap_bam
        Map[String, String]? rgid_2_fastq = MapRgid2Fqs.output_map

        Map[String, String]? runid_2_ont_basecall_model = GetBasecallModel.runid_2_model

        String last_postprocessing_date = today.yyyy_mm_dd
    }
}

task FixParticularBamHeaderIssue {
    meta {
        description: "Some one submitted to us BAMs with irregular headers. Fix that here."
    }
    input {
        File bam
    }
    output {
        File fixed = "~{prefix}.rg.fixed.bam"
    }

    Int disk_size = 100 + 2 * ceil(size(bam, 'GiB'))

    String prefix = basename(bam, '.bam')

    command <<<
        set -eux

        samtools view -H ~{bam} > original.header.txt
        cat original.header.txt

        # strip away the unwanted @RG line
        grep -vE "^@RG[[:space:]]SM:" original.header.txt > to.add.sample.name.header.txt
        diff original.header.txt to.add.sample.name.header.txt || true

        # add back the sample name to the correct @RG lines
        formatted_sample_name=$(grep -E "^@RG[[:space:]]SM:" original.header.txt | tr '\t' '\n' | grep "^SM:")
        TAB=$'\t'  # sed doesn't officially recoganize \t as tab
        for line_num in `grep -n "^@RG" to.add.sample.name.header.txt | awk -F ':' '{print $1}'`
        do
            echo "${line_num}"
            sed -i.bak "${line_num}s/$/""${TAB}""${formatted_sample_name}""/" to.add.sample.name.header.txt
        done
        cat to.add.sample.name.header.txt
        mv to.add.sample.name.header.txt fixed.header.txt

        samtools reheader fixed.header.txt ~{bam} > ~{prefix}.rg.fixed.bam
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk ~{disk_size} LOCAL"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
}