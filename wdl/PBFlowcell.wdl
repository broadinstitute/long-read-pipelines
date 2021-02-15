version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell.
## The workflow shards the subreads into clusters and performs CCS in parallel on each cluster.
## Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM.
## Various metrics are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow PBFlowcell {
    input {
        File bam
        File pbi

        String participant_name
        Int num_shards = 300
    }

    parameter_meta {
        bam:              "GCS path to raw subread bam"
        pbi:              "GCS path to pbi index for raw subread bam"

        participant_name: "name of the participant from whom these samples were obtained"
        num_shards:       "[default-valued] number of sharded BAMs to create (tune for performance)"
    }

    call PB.GetRunInfo { input: bam = bam }
    String ID = GetRunInfo.run_info["PU"]

    # break one raw BAM into fixed number of shards
    call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

    # then perform correction on each of the shard
    scatter (subreads in ShardLongReads.unmapped_shards) {
        call PB.CCS { input: subreads = subreads }
    }

    # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
    call Utils.MergeBams as MergeCorrected { input: bams = CCS.consensus, prefix = "~{participant_name}.~{ID}.corrected" }
    call PB.PBIndex as IndexCorrected { input: bam = MergeCorrected.merged_bam }
    call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report, prefix = "~{participant_name}.~{ID}" }

    File ccs_bam = MergeCorrected.merged_bam
    File ccs_pbi = IndexCorrected.pbi
    File ccs_report = MergeCCSReports.report

    call SummarizeCCSReport { input: report = ccs_report }

    output {
        File corrected_bam = ccs_bam
        File corrected_pbi = ccs_pbi
        File corrected_report = ccs_report

        Float zmws_input = SummarizeCCSReport.zmws_input
        Float zmws_pass_filters = SummarizeCCSReport.zmws_pass_filters
        Float zmws_fail_filters = SummarizeCCSReport.zmws_fail_filters
        Float zmws_pass_filters_pct = SummarizeCCSReport.zmws_pass_filters_pct
        Float zmws_fail_filters_pct = SummarizeCCSReport.zmws_fail_filters_pct
    }
}

task SummarizeCCSReport {
    input {
        File report

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(report, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{report} | grep 'ZMWs input' | awk -F": " '{ print $2 }' > zmws_input.txt
        cat ~{report} | grep 'ZMWs pass filters' | awk -F": " '{ print $2 }' | awk '{ print $1 }' > zmws_pass_filters.txt
        cat ~{report} | grep 'ZMWs fail filters' | awk -F": " '{ print $2 }' | awk '{ print $1 }' > zmws_fail_filters.txt
        cat ~{report} | grep 'ZMWs pass filters' | awk -F": " '{ print $2 }' | awk '{ print $2 }' | sed 's/[()%]//g' > zmws_pass_filters_pct.txt
        cat ~{report} | grep 'ZMWs fail filters' | awk -F": " '{ print $2 }' | awk '{ print $2 }' | sed 's/[()%]//g' > zmws_fail_filters_pct.txt
    >>>

    output {
        Float zmws_input = read_float("zmws_input.txt")
        Float zmws_pass_filters = read_float("zmws_pass_filters.txt")
        Float zmws_fail_filters = read_float("zmws_fail_filters.txt")
        Float zmws_pass_filters_pct = read_float("zmws_pass_filters_pct.txt")
        Float zmws_fail_filters_pct = read_float("zmws_fail_filters_pct.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.7"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
