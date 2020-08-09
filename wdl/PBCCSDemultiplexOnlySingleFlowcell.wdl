version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/ShardUtils.wdl" as SU
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/Finalize.wdl" as FF

workflow PBCCSDemultiplexOnlySingleFlowcell {
    input {
        String gcs_input_dir
        String? sample_name
        File barcode_file

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PU  = GetRunInfo.run_info["PU"]
        String ID  = PU
        String DIR = SM + "." + ID

        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = 2000000 }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
        }

        call Utils.MergeBams as MergeChunks { input: bams = CCS.consensus }
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }
    }

    if (length(FindBams.subread_bams) > 1) {
        call Utils.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }
    }

    File ccs_bam = select_first([ MergeRuns.merged_bam, MergeChunks.merged_bam[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    call PB.Demultiplex { input: bam = ccs_bam, prefix = "~{SM[0]}.~{ID[0]}", barcode_file = barcode_file }

    call PB.MakeSummarizedDemultiplexingReport as SummarizedDemuxReportPNG { input: report = Demultiplex.report }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeDemuxReads {
        input:
            files = Demultiplex.demux_bams,
            outdir = outdir + "/" + DIR[0] + "/reads"
    }

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ ccs_report ],
            outdir = outdir + "/" + DIR[0] + "/metrics/ccs"
    }

    call FF.FinalizeToDir as FinalizeLimaMetrics {
        input:
            files = [ Demultiplex.counts, Demultiplex.guess, Demultiplex.report, Demultiplex.summary ],
            outdir = outdir + "/" + DIR[0] + "/metrics/lima"
    }

    call FF.FinalizeToDir as FinalizeLimaSummary {
        input:
            files = SummarizedDemuxReportPNG.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/summary/png"
    }

    call FF.FinalizeToDir as FinalizeLimaDetailedPNG {
        input:
            files = DetailedDemuxReportPNG.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/detailed/png"
    }

    call FF.FinalizeToDir as FinalizeLimaDetailedPDF {
        input:
            files = DetailedDemuxReportPDF.report_files,
            outdir = outdir + "/" + DIR[0] + "/figures/lima/detailed/pdf"
    }
}

task CCS {
    input {
        String input_bam
        File input_bri
        File read_name_manifest

        Int min_passes = 3
        Float min_snr = 2.5
        Int min_length = 10
        Int max_length = 50000
        Float min_rq = 0.99
        Boolean by_strand = false

        Int batch_size = 10000

        Int cpus = 2
        Int mem = 16

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 4*ceil(size(input_bri, "GB"))

    command <<<
        set -ux

        # We'll batch our fetches into num_reads/batch_size requests.
        sed 's/,/\n/g' ~{read_name_manifest} | sort -t'/' -n -k2 > readnames.txt
        split -a 5 -d --additional-suffix=".txt" -l ~{batch_size} readnames.txt subchunk_

        # Get renewable auth token; see comment from @daviesrob at https://github.com/samtools/htslib/issues/803 .
        mkfifo /tmp/token_fifo
        ( while true ; do curl --retry 3 -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &

        # Fetch read batches in parallel, staggering requests a bit so that the Google access token website doesn't get
        # hammered with requests.  We'll also abort and restart if fetching a small batch of reads takes too long.
        HTS_AUTH_LOCATION=/tmp/token_fifo samtools view --no-PG -H ~{input_bam} > header.sam
        parallel --delay 30 --retries 3 --timeout 1800 'test $(cat {} | wc -l) -le $(cat {} | HTS_AUTH_LOCATION=/tmp/token_fifo bri get -i ~{input_bri} ~{input_bam} - | samtools view --no-PG -b > {.}.bam && samtools view {.}.bam | wc -l)' ::: subchunk_*

        # Merge all the pieces together.
        samtools cat --no-PG -h header.sam -o reads.bam subchunk_*.bam

        # Because we may process aligned data as well (where a read can appear more than once), we determine success
        # to be when the actual number of reads we extracted is equal to or greater than what we expected, rather
        # than strictly requiring each to be equal to one another.
        EXP_NUM_READS=$(cat readnames.txt | wc -l)
        ACT_NUM_READS=$(samtools view reads.bam | wc -l)
        echo $EXP_NUM_READS $ACT_NUM_READS

        if [ "$EXP_NUM_READS" -gt "$ACT_NUM_READS" ]
        then
            rm reads.bam
            exit 1
        fi

        ccs --min-passes ~{min_passes} \
            --min-snr ~{min_snr} \
            --min-length ~{min_length} \
            --max-length ~{max_length} \
            --min-rq ~{min_rq} \
            --num-threads ~{cpus} \
            --report-file ccs_report.txt \
            ~{if by_strand then "--by-strand" else ""} reads.bam ccs_unmapped.bam

        exit 0
    >>>

    output {
        File consensus = "ccs_unmapped.bam"
        File report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  7,
        max_retries:        3,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.22"
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

