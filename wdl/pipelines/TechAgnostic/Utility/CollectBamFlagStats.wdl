version 1.0

import "../../../tasks/Utility/BAMutils.wdl" as BU

workflow CollectBamFlagStats {
    meta {
        description: "Collect SAM flag stats of an aligned BAM"
    }
    input {
        File bam
        String disk_type = "HDD"
    }

    call BU.SamtoolsFlagStats  { input: bam = bam, output_format = 'JSON', disk_type = disk_type }
    call BU.ParseFlagStatsJson { input: sam_flag_stats_json = SamtoolsFlagStats.flag_stats }

    output {
        Map[String, Float] sam_flag_stats = ParseFlagStatsJson.qc_pass_reads_SAM_flag_stats
    }
}
