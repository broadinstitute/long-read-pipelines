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
        #File xml

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
    #call SummarizeXMLMetadata { input: xml = xml }

    call SummarizePBI as SummarizeSubreadsPBI { input: pbi = pbi }
    call SummarizePBI as SummarizeCCSPBI { input: pbi = ccs_pbi }
    call SummarizePBI as SummarizeCCSQ20PBI { input: pbi = ccs_pbi, qual_threshold = 20 }

    output {
        File corrected_bam = ccs_bam
        File corrected_pbi = ccs_pbi
        File corrected_report = ccs_report

        #Float xml_num_records = SummarizeXMLMetadata.xml_num_records
        #Float xml_total_length = SummarizeXMLMetadata.xml_total_length

        Float num_records = SummarizeSubreadsPBI.results['reads']
        Float total_length = SummarizeSubreadsPBI.results['bases']
        Float raw_yield = SummarizeSubreadsPBI.results['yield']

        Float polymerase_mean = SummarizeSubreadsPBI.results['polymerase_mean']
        Float polymerase_n50 = SummarizeSubreadsPBI.results['polymerase_n50']

        Float subread_mean = SummarizeSubreadsPBI.results['subread_mean']
        Float subread_n50 = SummarizeSubreadsPBI.results['subread_n50']

        Float ccs_num_records = SummarizeCCSPBI.results['reads']
        Float ccs_total_length = SummarizeCCSPBI.results['bases']
        Float ccs_mean_qual = SummarizeCCSPBI.results['mean_qual']
        Float ccs_yield = SummarizeCCSPBI.results['yield']

        Float ccs_num_records_q20 = SummarizeCCSQ20PBI.results['reads']
        Float ccs_total_length_q20 = SummarizeCCSQ20PBI.results['bases']
        Float ccs_mean_qual_q20 = SummarizeCCSQ20PBI.results['mean_qual']
        Float ccs_yield_q20 = SummarizeCCSQ20PBI.results['yield']

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

task SummarizeXMLMetadata {
    input {
        File xml

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(xml, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{xml} | grep '<pbds:TotalLength>' | sed 's/<pbds:TotalLength>//g' | sed 's/<\/pbds:TotalLength>//' | sed 's/\s*//g' > xml_total_length.txt
        cat ~{xml} | grep '<pbds:NumRecords>' | sed 's/<pbds:NumRecords>//g' | sed 's/<\/pbds:NumRecords>//' | sed 's/\s*//g' > xml_num_records.txt
    >>>

    output {
        Float xml_total_length = read_float("xml_total_length.txt")
        Float xml_num_records = read_float("xml_num_records.txt")
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

task SummarizePBI {
    input {
        File pbi
        Int qual_threshold = 0

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(pbi, "GB"))

    command <<<
        set -euxo pipefail

        python3 /usr/local/bin/compute_pbi_stats.py -q ~{qual_threshold} ~{pbi} > map.txt

        cat map.txt
    >>>

    output {
        Map[String, Int] results = read_map("map.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.25"
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
