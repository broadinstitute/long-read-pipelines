version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/ShardUtils.wdl" as SU
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Finalize.wdl" as FF
import "tasks/CallSVs.wdl" as SV
import "tasks/CallSmallVariants.wdl" as SMV

workflow PBCCSDemultiplexWholeGenomeSingleFlowcell {
    input {
        String gcs_input_dir
        String? sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File tandem_repeat_bed
        File ref_flat
        File dbsnp_vcf
        File dbsnp_tbi

        String mt_chr_name
        File metrics_locus

        Int num_shards = 300

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID
        String RG  = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call SU.IndexUnalignedBam { input: input_bam = subread_bam }
        call SU.MakeReadNameManifests { input: input_bri = IndexUnalignedBam.bri, N = num_shards }

        scatter (manifest in MakeReadNameManifests.manifest_chunks) {
            call CCS {
                input:
                    input_bam          = subread_bam,
                    input_bri          = IndexUnalignedBam.bri,
                    read_name_manifest = manifest,
            }
        }

        call AR.MergeBams as MergeChunks { input: bams = CCS.consensus }
        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }
    }

    if (length(FindBams.subread_bams) > 1) {
        call AR.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }
    }

    File ccs_bam = select_first([ MergeRuns.merged_bam, MergeChunks.merged_bam[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    call PB.Demultiplex { input: bam = ccs_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    call PB.MakeSummarizedDemultiplexingReport as SummarizedDemuxReportPNG { input: report = Demultiplex.report }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPNG { input: report = Demultiplex.report, type="png" }
    call PB.MakeDetailedDemultiplexingReport as DetailedDemuxReportPDF { input: report = Demultiplex.report, type="pdf" }

    scatter (demux_bam in Demultiplex.demux_bams) {
        String BC   = sub(basename(demux_bam, ".bam"), "~{SM[0]}.~{ID[0]}.", "")
        String BCID = ID[0] + "." + BC
        String BCRG = "@RG\\tID:~{BCID}\\tSM:~{SM[0]}\\tPL:~{PL[0]}\\tPU:~{PU[0]}\\tDT:~{DT[0]}"

        call AR.Minimap2 as AlignBarcode {
            input:
                reads      = [ demux_bam ],
                ref_fasta  = ref_fasta,
                RG         = BCRG,
                map_preset = "asm20",
                prefix     = BCID + ".aligned"
        }

        call AM.AlignedMetrics as PerBarcodeMetrics {
            input:
                aligned_bam    = AlignBarcode.aligned_bam,
                aligned_bai    = AlignBarcode.aligned_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                ref_flat       = ref_flat,
                dbsnp_vcf      = dbsnp_vcf,
                dbsnp_tbi      = dbsnp_tbi,
                metrics_locus  = metrics_locus,
                per            = "flowcell",
                type           = "barcode",
                label          = BC,
                gcs_output_dir = outdir + "/" + DIR[0]
        }

        call Utils.BamToBed { input: bam = AlignBarcode.aligned_bam, prefix = "~{SM[0]}.~{ID[0]}.~{BC}" }

        call SV.CallSVs as CallSVs {
            input:
                bam               = AlignBarcode.aligned_bam,
                bai               = AlignBarcode.aligned_bai,

                ref_fasta         = ref_fasta,
                ref_fasta_fai     = ref_fasta_fai,
                tandem_repeat_bed = tandem_repeat_bed
        }

        call SMV.CallSmallVariants as CallSmallVariants {
            input:
                bam               = AlignBarcode.aligned_bam,
                bai               = AlignBarcode.aligned_bai,

                ref_fasta         = ref_fasta,
                ref_fasta_fai     = ref_fasta_fai,
                ref_dict          = ref_dict
        }

        call FF.FinalizeToDir as FinalizeDemuxAlignedReads {
            input:
                files  = [ AlignBarcode.aligned_bam, AlignBarcode.aligned_bai, BamToBed.bed ],
                outdir = outdir + "/" + DIR[0] + "/alignments/" + BC
        }

        call FF.FinalizeToDir as FinalizeSVs {
            input:
                files = [ CallSVs.sniffles_vcf, CallSVs.svim_vcf ],
                outdir = outdir + "/" + DIR[0] + "/variants/" + BC
        }

        call FF.FinalizeToDir as FinalizeSmallVariants {
            input:
                files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
                outdir = outdir + "/" + DIR[0] + "/variants/" + BC
        }
    }

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
            reads.bam \
            ccs_unmapped.bam

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

