version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell.
## The workflow shards the subreads into clusters and performs CCS in parallel on each cluster.
## Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM.
## Various metrics are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/NanoPlot.wdl" as NP
import "tasks/Finalize.wdl" as FF
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/SMRTtools.wdl"
import "tasks/metrics/CollectSMRTCellUnalignedMetrics.wdl" as uBAMCustomMetrics

workflow PBNonBarcodedWGSSMRTcellPreprocess {
    input {
        String smrtcell_data_dir
        String experiment_type

        String SM
        String LB

        Int num_shards = 50

        String gcs_out_root_dir

        Boolean validate_shards = false
    }

    parameter_meta {
        smrtcell_data_dir:  "Bucket location of the 'folder' that contains the CCSed BAM, PBI, XMLs and other companion files"

        SM:                 "the value to place in the BAM read group's SM field"
        LB:                 "the value to place in the BAM read group's LB (library) field"

        num_shards:         "[default-valued] number of shards into which fastq files should be batched"
        experiment_type:    "type of experiment run (CLR, CCS, ISOSEQ, MASSEQ)"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        validate_shards: "an internal control step to validate data coherency; turning this on increases wait time"
    }

    String workflow_name = "PBNonBarcodedWGSSMRTcellPreprocess"

    # prep
    call DistinguishFiles as prep { input: smrtcell_data_dir = smrtcell_data_dir }
    call PB.GetRunInfo { input: bam = bam, SM = SM }
    String PU = GetRunInfo.run_info['PU']
    String bam = prep.bam_path
    String pbi = prep.pbi_path
    String movie = basename(basename(basename(prep.bam_path, ".reads.bam"), ".hifi_reads.bam"), ".subreads.bam")

    # runqc-reports
    call SMRTtools.RunQCreports { input: smrtcell_data_dir = smrtcell_data_dir, movie_name = movie }

    # collect metrics on raw PBI
    call uBAMCustomMetrics.CollectSMRTCellUnalignedMetrics {input: smrtcell_pbi = prep.pbi_path}

    # then perform correction on each of the shard, if applicable and data hasn't been
    call Utils.ComputeAllowedLocalSSD as Guess {input: intended_gb = 3*ceil(size(bam, "GB") + size(pbi, "GB"))}
    call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}
    call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards, num_ssds = Guess.numb_of_local_ssd, zones = arbitrary.zones}
    scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
        # sometimes we see the sharded bams mising EOF marker, use this as
        if (validate_shards) {call Utils.CountBamRecords as ValidateShard {input: bam = unmapped_shard}}

        if (experiment_type != "CLR") {
            if (!GetRunInfo.is_corrected) { call PB.CCS { input: subreads = unmapped_shard } }
            call PB.ExtractHifiReads {
                input:
                    bam = select_first([CCS.consensus, unmapped_shard]),
                    sample_name = SM,
                    library     = LB,
                    prefix      = "~{movie}.hifi_reads"
            }
        }

        File unaligned_bam = select_first([ExtractHifiReads.hifi_bam, unmapped_shard])
        call Utils.BamToFastq { input: bam = unaligned_bam, prefix = basename(unaligned_bam, ".bam") }
    }
    # FASTQs: hifi_fastq if ccs, or raw fastq otherwise
    call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = BamToFastq.reads_fq }

    # if applicable: merge corrected, unaligned reads, and pipeline-CCSed reports
    if (experiment_type != "CLR") {
        call Utils.MergeBams as MergeHiFiUnalignedReads { input: bams = select_all(ExtractHifiReads.hifi_bam), prefix = "~{PU}.reads" }
        call PB.PBIndex as IndexHiFiUnalignedReads { input: bam = MergeHiFiUnalignedReads.merged_bam }

        if (!GetRunInfo.is_corrected) {
            call PB.MergeCCSReports as MergeCCSReports { input: reports = select_all(CCS.report), prefix = PU }
        }

        call PB.SummarizeCCSReport { input: report = select_first([prep.ccs_report_txt_path, MergeCCSReports.report]) }
    }

    ###################################################################################
    # Finalize data
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/" + movie
    String outdir_metrics = outdir + "/metrics"
    # runqc-reports
    call FF.FinalizeToFile as FinalizeRunQCreports {
        input:
            file = RunQCreports.runqc_reports_tar_gz,
            outdir = outdir_metrics
    }
    # metrics on raw PBI
    call FF.FinalizeToFile as FinalizePBISummary {
        input:
            file = CollectSMRTCellUnalignedMetrics.pbi_summary,
            outdir = outdir_metrics
    }
    # FASTQ
    call FF.FinalizeToFile as FinalizeFastq {
        input:
            outdir = outdir,
            file = MergeAllFastqs.merged_fastq,
            name = movie + if (experiment_type != "CLR") then ".hifi.fq.gz" else ".subreads.fq.gz"
    }
    # if applicable: hifi bam/pbi, ccs report
    if (experiment_type != "CLR") {
        call FF.FinalizeToFile as FinalizeHiFiUnalignedBam { input: outdir = outdir, file = select_first([MergeHiFiUnalignedReads.merged_bam]) }
        call FF.FinalizeToFile as FinalizeHiFiUnalignedPbi { input: outdir = outdir, file = select_first([IndexHiFiUnalignedReads.pbi]) }

        if (!GetRunInfo.is_corrected) {
            call FF.FinalizeToFile as FinalizeCCSReport { input: outdir = outdir, file = select_first([MergeCCSReports.report]) }
        }
    }

    call GU.GetTodayDate as today {}

    output {
        # files for downstream WDLs
        File fastq = FinalizeFastq.gcs_path # not necessarily hifi
        File? hifi_bam = FinalizeHiFiUnalignedBam.gcs_path
        File? hifi_pbi = FinalizeHiFiUnalignedPbi.gcs_path
        File? ccs_report = FinalizeCCSReport.gcs_path

        # metrics
        File runqc_reports_tar_gz = FinalizeRunQCreports.gcs_path

        File smrtcell_pbi_summary = FinalizePBISummary.gcs_path

        Float? zmws_ccs_input = SummarizeCCSReport.zmws_input
        Float? zmws_ccs_pass_filters = SummarizeCCSReport.zmws_pass_filters
        Float? zmws_ccs_fail_filters = SummarizeCCSReport.zmws_fail_filters
        Float? zmws_ccs_shortcut_filters = SummarizeCCSReport.zmws_shortcut_filters
        Float? zmws_ccs_pass_filters_pct = SummarizeCCSReport.zmws_pass_filters_pct
        Float? zmws_ccs_fail_filters_pct = SummarizeCCSReport.zmws_fail_filters_pct
        Float? zmws_ccs_shortcut_filters_pct = SummarizeCCSReport.zmws_shortcut_filters_pct

        String last_preprocessing_date = today.yyyy_mm_dd
    }
}

task DistinguishFiles {
    meta {
        desciption: "Given a path to an CCSed SMRTCell's data folder mirrored onto the bucket, separate files for different downstream processing"
    }

    input {
        String smrtcell_data_dir
    }

    command <<<
        set -eux

        gsutil ls ~{smrtcell_data_dir} > files.list

        grep -E "consensusreadset.xml$" files.list > main.xml.txt
        grep -E "sts.xml$" files.list > sts.xml.txt
        grep -E "reads.bam$" files.list > bam.txt
        grep -E ".pbi$" files.list > pbi.txt

        grep -E "ccs_reports.txt$" files.list > ccs_report.txt

        if grep -qF 'hifi' bam.txt; then echo "true" > is_hifi.txt; else echo "false" > is_hifi.txt; fi
    >>>

    output {
        String main_xml_path = read_string("main.xml.txt")
        String bam_path = read_string("bam.txt")
        String pbi_path = read_string("pbi.txt")
        String? ccs_report_txt_path = read_string("ccs_report.txt")

        Boolean is_hifi = read_boolean("is_hifi.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
