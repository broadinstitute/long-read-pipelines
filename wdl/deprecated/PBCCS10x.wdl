version 1.0

import "../tasks/Utility/PBUtils.wdl" as PB
import "../tasks/Utility/Utils.wdl" as Utils
import "../tasks/Alignment/AlignReads.wdl" as AR
import "../tasks/QC/AlignedMetrics.wdl" as AM
import "tasks/AnnotateAdapters.wdl" as AA
import "../tasks/Utility/Finalize.wdl" as FF

workflow PBCCS10x {
    input {
        Array[File] bams
        File ref_map_file

        String participant_name
        File barcode_file
        Int num_shards = 300

        Boolean drop_per_base_N_pulse_tags = true

        String gcs_out_root_dir
    }

    parameter_meta {
        bams:             "GCS path to raw subreads or CCS data"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"

        participant_name: "name of the participant from whom these samples were obtained"
        barcode_file:     "GCS path to the fasta file that specifies the expected set of multiplexing barcodes"
        num_shards:       "[default-valued] number of sharded BAMs to create (tune for performance)"

        gcs_out_root_dir: "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCS10x/" + participant_name

    # scatter over all sample BAMs
    scatter (bam in bams) {
        File pbi = sub(bam, ".bam$", ".bam.pbi")

        call PB.GetRunInfo { input: bam = bam }
        String ID = GetRunInfo.run_info["PU"]

        # break one raw BAM into fixed number of shards
        call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

        # then perform correction on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }

            call AA.AnnotateAdapters { input: bam = CCS.consensus }

            call PB.Align as AlignCorrected {
                input:
                    bam         = AnnotateAdapters.annotated_bam,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = participant_name,
                    drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags,
                    map_preset  = "ISOSEQ"
            }

            call Utils.CountBamRecords as CountSubreadsInShard { input: bam = subreads }
            call Utils.CountBamRecords as CountCorrectedReadsInShard { input: bam = CCS.consensus }
            call Utils.CountBamRecords as CountAnnotatedReadsInShard { input: bam = AnnotateAdapters.annotated_bam }
        }

#        call C3.Cat as CountNumPassesInRun { input: files = CountNumPasses.num_passes, out = "num_passes.txt" }

        call Utils.Sum as CountSubreadsInFlowcell { input: ints = CountSubreadsInShard.num_records }
        call Utils.Sum as CountAnnotatedReadsInFlowcell { input: ints = CountAnnotatedReadsInShard.num_records }
        call Utils.Sum as CountCorrectedReadsInFlowcell { input: ints = CountCorrectedReadsInShard.num_records }

        call Utils.MergeBams as MergeAnnotated { input: bams = AnnotateAdapters.annotated_bam }
        call Utils.MergeBams as MergeCorrected { input: bams = AlignCorrected.aligned_bam }

        # compute alignment metrics
        call AM.AlignedMetrics as PerFlowcellMetrics {
            input:
                aligned_bam    = MergeCorrected.merged_bam,
                aligned_bai    = MergeCorrected.merged_bai,
                ref_fasta      = ref_map['fasta'],
                ref_dict       = ref_map['dict'],
                gcs_output_dir = outdir + "/metrics/per_flowcell/" + ID
        }

        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report, prefix = "~{participant_name}.~{ID}" }

        call FF.FinalizeToDir as FinalizeCCSReport {
            input:
                files = [ MergeCCSReports.report ],
                outdir = outdir + "/metrics/per_flowcell/" + ID + "/ccs_metrics"
        }
    }

#    call C3.Cat as CountNumPassesAll { input: files = CountNumPassesInRun.merged, out = "num_passes.txt" }

    call Utils.Sum as CountSubreads { input: ints = CountSubreadsInFlowcell.sum, prefix = "num_subreads" }
    call Utils.Sum as CountCorrectedReads { input: ints = CountCorrectedReadsInFlowcell.sum, prefix = "num_consensus" }
    call Utils.Sum as CountAnnotatedReads { input: ints = CountAnnotatedReadsInFlowcell.sum, prefix = "num_annotated" }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllCorrected { input: bams = MergeCorrected.merged_bam, prefix = "~{participant_name}.corrected" }
        call Utils.MergeBams as MergeAllAnnotated { input: bams = MergeAnnotated.merged_bam, prefix = "~{participant_name}.annotated" }

        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }
    }

    File ccs_bam = select_first([ MergeAllCorrected.merged_bam, MergeCorrected.merged_bam[0] ])
    File ccs_bai = select_first([ MergeAllCorrected.merged_bai, MergeCorrected.merged_bai[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    File annotated_bam = select_first([ MergeAllAnnotated.merged_bam, MergeAnnotated.merged_bam[0] ])
    File annotated_bai = select_first([ MergeAllAnnotated.merged_bai, MergeAnnotated.merged_bai[0] ])

    call GrepCountBamRecords as GrepAnnotatedReadsWithCBC {
        input:
            bam = annotated_bam,
            prefix = "num_annotated_with_cbc",
            regex = "CB:Z:[ACGT]"
    }

    call GrepCountBamRecords as GrepAnnotatedReadsWithCBCAndUniqueAlignment {
        input:
            bam = annotated_bam,
            samfilter = "-F 0x100",
            prefix = "num_annotated_with_cbc_and_unique_alignment",
            regex = "CB:Z:[ACGT]"
    }

    call BamToTable { input: bam = annotated_bam, prefix = "reads_aligned_annotated.table" }

    # compute alignment metrics
    call AM.AlignedMetrics as PerSampleMetrics {
        input:
            aligned_bam    = annotated_bam,
            aligned_bai    = annotated_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeCorrectedReadCounts {
        input:
            files = [ CountSubreads.sum_file, CountAnnotatedReads.sum_file, CountCorrectedReads.sum_file,
                      GrepAnnotatedReadsWithCBC.num_records_file, GrepAnnotatedReadsWithCBCAndUniqueAlignment.num_records_file ],
            outdir = outdir + "/metrics/read_counts"
    }

#    call FF.FinalizeToDir as FinalizeNumPasses {
#        input:
#            files = [ CountNumPassesAll.merged ],
#            outdir = outdir + "/metrics/num_passes"
#    }

    call FF.FinalizeToDir as FinalizeBamTable {
        input:
            files = [ BamToTable.table ],
            outdir = outdir + "/metrics/bam_tables"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ annotated_bam, annotated_bai ],
            outdir = outdir + "/alignments"
    }
}

task GrepCountBamRecords {

    meta {
        description : "Count the number of records in a bam file that match a given regex."
    }

    parameter_meta {
        bam: "BAM file to be filtered."
        samfilter: "[Optional] Extra arguments to pass into samtools view."
        regex: "Regex to match against the bam file."
        invert: "[Optional] Invert the regex match."
        prefix: "[Optional] Prefix string to name the output file (Default: sum)."
    }

    input {
        File bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-vc" else "-c"

    command <<<
        set -euxo pipefail

        samtools view ~{samfilter} ~{bam} | grep ~{arg} ~{regex} > ~{prefix}.txt
    >>>

    output {
        Int num_records = read_int("~{prefix}.txt")
        File num_records_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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

task BamToTable {

    meta {
        description : "Convert a BAM file to a table txt file."
    }

    parameter_meta {
        bam: "BAM file to be converted."
        prefix: "Prefix for the output table."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam, "GB"))

    command <<<
        samtools view ~{bam} | perl -n -e '($nm) = $_ =~ /NM:i:(\d+)/; ($as) = $_ =~ /AS:i:(\d+)/; ($za) = $_ =~ /ZA:Z:(\w+|\.)/; ($zu) = $_ =~ /ZU:Z:(\w+|\.)/; ($cr) = $_ =~ /CR:Z:(\w+|\.)/; ($cb) = $_ =~ /CB:Z:(\w+|\.)/; @a = split(/\s+/); print join("\t", $a[0], $a[1], $a[2], $a[3], $a[4], length($a[9]), $nm, $as, $za, $zu, $cr, $cb, $a[1], ($a[1] & 0x1 ? "paired" : "unpaired"), ($a[1] & 0x4 ? "unmapped" : "mapped"), ($a[1] & 0x10 ? "rev" : "fwd"), ($a[1] & 0x100 ? "secondary" : "primary"), ($a[1] & 0x800 ? "supplementary" : "non_supplementary")) . "\n"' | gzip > ~{prefix}.txt.gz
    >>>

    output {
        File table = "~{prefix}.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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
