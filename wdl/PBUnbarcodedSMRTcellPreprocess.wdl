version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/NanoPlot.wdl"

import "tasks/Finalize.wdl" as FF

workflow PBUnbarcodedSMRTcellPreprocess {
    meta {
        description: "Performs CCS on old SMRT cells that were not on-board CCS processed."
    }
    input {
        File? subreads_bam
        File? subreads_pbi
        String? smrt_folder

        String movie

        String sample
        String participant_name

        Int? num_shards

        String gcs_out_root_dir
    }

    String workflow_name = "PBUnbarcodedSMRTcellPreprocess"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{participant_name}/~{workflow_name}"
    String output_prefix = sample + "." + movie

    ###################################################################################
    # check inputs
    Boolean incompatible_inputs = false
    if (defined(smrt_folder)) {
        if (defined(subreads_bam) || defined(subreads_pbi)) {
            call Utils.StopWorkflow as IncompatibleInputs {input: reason = "Please specify only smrt_folder, or (subreads_bam, subreads_pbi), not both."}
        }
    }
    if (!defined(smrt_folder)) {
        if (!(defined(subreads_bam) && defined(subreads_pbi))) {
            call Utils.StopWorkflow as PartialInput {input: reason = "Please specify both (subreads_bam, subreads_pbi), when input is subreads of CCS library."}
        }
    }

    ###################################################################################
    Boolean need_cloud_ccs = !defined(smrt_folder)
    if (need_cloud_ccs) {
        call PB.SummarizePBI as SummarizeSubreadsPBI { input: pbi = select_first([subreads_pbi]), runtime_attr_override = { 'mem_gb': 72 } }
    }
    if (!need_cloud_ccs) {
        call DistinguishFiles {input: smrtcell_data_dir = select_first([smrt_folder])}
    }

    File raw_bam = select_first([subreads_bam, DistinguishFiles.bam_path])
    File raw_pbi = select_first([subreads_pbi, DistinguishFiles.pbi_path])

    if (need_cloud_ccs) {
        ##########
        # shard-[ccs,...]-merge
        ##########
        call Utils.ComputeAllowedLocalSSD as Guess {input: intended_gb = ceil(size(raw_bam, "GB") * 4 + size(raw_pbi, "GB"))}
        call Utils.RandomZoneSpewer as arbitrary {input: num_of_zones = 3}
        call PB.ShardLongReads as ShardLongReads {
            input:
                unaligned_bam = raw_bam,
                unaligned_pbi = raw_pbi,
                num_shards = num_shards,
                drop_per_base_N_pulse_tags = true,
                num_ssds = Guess.numb_of_local_ssd,
                zones = arbitrary.zones
        }

        # then perform correction and alignment on each of the shard
        scatter (shard in ShardLongReads.unmapped_shards) {
            # sometimes we see the sharded bams mising EOF marker, use this as
            Boolean validate_shards = true
            if (validate_shards) {call Utils.CountBamRecords as ValidateShard {input: bam = shard}}

            call PB.CCS as CCS { input: subreads = shard }

            call PB.ExtractHifiReads as ExtractFromCloudCCSShard {
                input:
                    bam = CCS.consensus,
                    sample_name = sample,
                    library     = movie
            }
            call Utils.BamToFastq as BamToFastqCloudCCS { input: bam = ExtractFromCloudCCSShard.hifi_bam, prefix = "doesnotmatter" }
        }

        # Merge
        call Utils.MergeFastqs  as MergeAllFastqs         { input: fastqs  = BamToFastqCloudCCS.reads_fq,       prefix = output_prefix }
        call PB.MergeCCSReports as MergeCCSReports        { input: reports = select_all(CCS.report),            prefix = output_prefix }
        call Utils.MergeBams    as MergeCCSUnalignedReads { input: bams    = ExtractFromCloudCCSShard.hifi_bam, prefix = output_prefix }
    }

    if (!need_cloud_ccs) {
        call PB.ExtractHifiReads as ExtractFromOnBoardCCS {
            input:
                bam = select_first([DistinguishFiles.bam_path]),
                sample_name = sample,
                library     = movie
        }
        call Utils.BamToFastq as BamToFastqOnBoardCCS { input: bam = ExtractFromOnBoardCCS.hifi_bam, prefix = "doesnotmatter" }
    }

    File ccs_report = select_first([MergeCCSReports.report, DistinguishFiles.ccs_report_txt_path])
    File extracted_hifi_bam = select_first([MergeCCSUnalignedReads.merged_bam, ExtractFromOnBoardCCS.hifi_bam])
    File extracted_hifi_fastq = select_first([MergeAllFastqs.merged_fastq, BamToFastqOnBoardCCS.reads_fq])

    # no BAI because that's not quite useful here
    call PB.PBIndex as IndexHiFiReads { input: bam = extracted_hifi_bam }

    # quick metrics collection
    call PB.SummarizeCCSReport as SummarizeCCSReport   { input: report = ccs_report }
    call PB.SummarizePBI       as SummarizeHiFiReadsPBI { input: pbi = IndexHiFiReads.pbi, runtime_attr_override = { 'mem_gb': 72 } }
    call NanoPlot.NanoPlotFromUBam {input: bam = extracted_hifi_bam}

    ###################################################################################
    ##########
    # store the results into designated bucket
    ##########
    String cdir = outdir
    File keyfile = IndexHiFiReads.pbi

    if (need_cloud_ccs) {
        call FF.FinalizeToFile as FinalizeCCSReport {
            input: outdir = cdir, file = ccs_report, keyfile = keyfile
        }
    }

    call FF.FinalizeToFile as FinalizeHiFiFastq {
        input: outdir = cdir, file = extracted_hifi_fastq, name = output_prefix + ".hifi_reads.fastq.gz", keyfile = keyfile
    }

    call FF.FinalizeToFile as FinalizeHiFiBam {
        input: outdir = cdir, file = extracted_hifi_bam, name = output_prefix + ".hifi_reads.bam", keyfile = keyfile
    }
    call FF.FinalizeToFile as FinalizeHiFiPbi {
        input: outdir = cdir, file = IndexHiFiReads.pbi, name = output_prefix + ".hifi_reads.bam.pbi", keyfile = keyfile
    }

    call GU.GetTodayDate as today {}

    output {
        String last_preprocess_date = today.yyyy_mm_dd

        # hifi reads related
        File hifi_fastq = FinalizeHiFiFastq.gcs_path
        File hifi_bam   = FinalizeHiFiBam.gcs_path
        File hifi_pbi   = FinalizeHiFiPbi.gcs_path
        Map[String, Float] hifi_reads_stats = NanoPlotFromUBam.stats_map

        # ccs reads related
        Map[String, Float] ccs_reads_stats = SummarizeCCSReport.stats

        # subreads related, but is output only when cloud-ran CCS happened
        File? cloud_gen_ccs_report = FinalizeCCSReport.gcs_path
        Map[String, Float]? subreads_stats = SummarizeSubreadsPBI.results
    }
}

task DistinguishFiles {
    meta {
        desciption: "CAUTION: only works for SMRTLink v10, v11. Given a path to an CCSed SMRTCell's data folder mirrored onto the bucket, separate files for different downstream processing"
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
        String ccs_report_txt_path = read_string("ccs_report.txt")

        Boolean is_hifi = read_boolean("is_hifi.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
