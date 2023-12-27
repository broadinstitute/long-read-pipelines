version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/BAMutils.wdl" as BU

workflow CountTheBeans {
    meta {
        desciption: "For long-read bams (PacBio and ONT), gather information about reads with and without the ML/MM tags for methylation."
    }
    input {
        File  bam
        File? bai
        String bam_descriptor
        Boolean use_local_ssd = false
        String? gcs_out_root_dir
    }


    call BU.CountMethylCallReads as Count { input: bam = bam, bai = bai, disk_type = if(use_local_ssd) then "LOCAL" else "SSD"}

    call BU.GatherReadsWithoutMethylCalls as GatherBitter { input: bam = bam, bai = bai, disk_type = if(use_local_ssd) then "LOCAL" else "SSD"}

    if (defined(gcs_out_root_dir)) {
        String gcs_out = sub(select_first([gcs_out_root_dir]), "/$", "") + "/"
        String out_dir = gcs_out + basename(bam, ".bam") + "." + bam_descriptor + "/RecordsWithoutMethylTags/"
        call FF.FinalizeToDir {
            input:
                files = [GatherBitter.no_ml_reads,
                         GatherBitter.no_mm_reads,
                         GatherBitter.names_missing_only_one_tag,
                         GatherBitter.names_missing_both_tags],
                outdir = out_dir
        }
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
}
