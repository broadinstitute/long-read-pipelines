version 1.0

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallSVs.wdl" as SV
import "tasks/CallSmallVariants.wdl" as SMV
import "tasks/Figures.wdl" as FIG
import "tasks/Methylation.wdl" as Meth
import "tasks/Finalize.wdl" as FF

workflow ONTFlowcell {
    input {
        File final_summary
        File sequencing_summary
        File ref_map_file

        String participant_name
        Int num_shards = 50

        String gcs_out_root_dir
    }

    parameter_meta {
        final_summary:      "GCS path to '*final_summary*.txt*' file for basecalled fastq files"
        sequencing_summary: "GCS path to '*sequencing_summary*.txt*' file for basecalled fastq files"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        participant_name:   "name of the participant from whom these samples were obtained"
        num_shards:         "[default-valued] number of shards into which fastq files should be batched"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTFlowcell/" + participant_name

    call ONT.GetRunInfo { input: summary_file = final_summary }
    call ONT.ListFiles as ListFast5s { input: summary_file = final_summary, suffix = "fast5" }
    call ONT.ListFiles as ListFastqs { input: summary_file = final_summary, suffix = "fastq" }

    String SM  = participant_name
    String PL  = "ONT"
    String PU  = GetRunInfo.run_info["instrument"]
    String DT  = GetRunInfo.run_info["started"]
    String ID  = GetRunInfo.run_info["flow_cell_id"] + "." + GetRunInfo.run_info["position"]
    String DIR = GetRunInfo.run_info["protocol_group_id"] + "." + SM + "." + ID
    String SID = ID + "." + sub(GetRunInfo.run_info["protocol_run_id"], "-.*", "")
    String RG = "@RG\\tID:~{SID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

    call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = ListFastqs.manifest, N = num_shards }

    scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
        call AR.Minimap2 as AlignReads {
            input:
                reads      = read_lines(manifest_chunk),
                ref_fasta  = ref_map['fasta'],
                RG         = RG,
                map_preset = "map-ont"
        }
    }

    call Utils.MergeBams as MergeReads { input: bams = AlignReads.aligned_bam, prefix = "~{participant_name}.~{ID}" }

    call AM.AlignedMetrics as PerFlowcellMetrics {
        input:
            aligned_bam    = MergeReads.merged_bam,
            aligned_bai    = MergeReads.merged_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics/per_flowcell/" + SID
    }

    call FIG.Figures as PerFlowcellFigures {
        input:
            summary_files  = [ sequencing_summary ],
            gcs_output_dir = outdir + "/metrics/per_flowcell/" + SID
    }

    call SummarizeNanoStats { input: report = PerFlowcellFigures.NanoPlotFromSummaryStats }

    output {
        File aligned_bam = MergeReads.merged_bam
        File aligned_bai = MergeReads.merged_bai

        Float active_channels = SummarizeNanoStats.results['Active_channels']

        Float num_reads = SummarizeNanoStats.results['Number_of_reads']
        Float total_bases = SummarizeNanoStats.results['Total_bases']
        Float raw_yield = SummarizeNanoStats.results['Total_bases']/3088286401.0

        Float read_length_mean = SummarizeNanoStats.results['Mean_read_length']
        Float read_length_N50 = SummarizeNanoStats.results['Read_length_N50']
        Float read_length_median = SummarizeNanoStats.results['Median_read_length']

        Float read_qual_mean = SummarizeNanoStats.results['Mean_read_quality']
        Float read_qual_median = SummarizeNanoStats.results['Mean_read_quality']

        Float num_reads_gt_Q5 = SummarizeNanoStats.results['Number_of_reads_gt_Q5']
        Float num_reads_gt_Q7 = SummarizeNanoStats.results['Number_of_reads_gt_Q7']
        Float num_reads_gt_Q10 = SummarizeNanoStats.results['Number_of_reads_gt_Q10']
        Float num_reads_gt_Q12 = SummarizeNanoStats.results['Number_of_reads_gt_Q12']
        Float num_reads_gt_Q15 = SummarizeNanoStats.results['Number_of_reads_gt_Q15']
    }
}

task SummarizeNanoStats {
    input {
        File report

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(report, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{report} | grep -v -e '^General' -e '^Number,' | head -13 | sed 's/:\s\+/\t/g' | sed 's/,//g' | sed 's/ /_/g' | awk '{ print $1 "\t" $2 }' | sed 's/_(.*//g' | sed 's/>/Number_of_reads_gt_/' > map.txt
        cat map.txt
    >>>

    output {
        Map[String, Float] results = read_map("map.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-ont:0.1.1"
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
