version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/Utils.wdl"

workflow CountTheBeans {
    meta {
        desciption:
        "For long-read bams (PacBio and ONT), gather information about reads with and without the ML/MM tags for methylation."
    }

    parameter_meta {
        bam_descriptor:
        "a short description (no space) of the purpose of the BAM; used only for saving the results."
        save_read_names_only:
        "the workflow extracts both the reads and the read names missing the methylation SAM tags; if this is set to true, then only the read names are saved."
    }
    input {
        File  bam
        File? bai

        Boolean use_local_ssd = false

        String? bam_descriptor
        Boolean save_read_names_only = true
        String? gcs_out_root_dir
    }

    output {
        Map[String, String] methyl_tag_simple_stats = {
                                        'raw_record_cnt': Count.raw_count,
                                        'raw_record_with-mm-ml_cnt': Count.bean_count,
                                        'primary_record_cnt': Count.non_2304_count,
                                        'primary_record_with-mm-ml_cnt': Count.non_2304_bean_count,
                                        'files_holding_reads_without_tags': select_first([FinalizeToDir.gcs_dir, "None"])
        }
    }

    if (defined(gcs_out_root_dir)!=defined(bam_descriptor)) {
        call Utils.StopWorkflow { input: reason = "'bam_descriptor' and 'gcs_out_root_dir' must be provided/omitted together." }
    }

    call BU.CountMethylCallReads as Count { input: bam = bam, bai = bai, disk_type = if(use_local_ssd) then "LOCAL" else "SSD"}

    call BU.GatherReadsWithoutMethylCalls as GatherBitter { input: bam = bam, bai = bai, disk_type = if(use_local_ssd) then "LOCAL" else "SSD"}

    if (defined(gcs_out_root_dir)) {
        String gcs_out = sub(select_first([gcs_out_root_dir]), "/$", "") + "/"
        String out_dir = gcs_out + basename(bam, ".bam") + "." + select_first([bam_descriptor]) + "/RecordsWithoutMethylTags/"

        Array[File] reads = [GatherBitter.no_ml_reads, GatherBitter.no_mm_reads]
        Array[File] read_names = [GatherBitter.names_missing_only_one_tag, GatherBitter.names_missing_both_tags]
        Array[File] files_to_save = if (save_read_names_only) then read_names else flatten([reads, read_names])

        call FF.FinalizeToDir { input:
            files = files_to_save,
            outdir = out_dir
        }
    }
}
