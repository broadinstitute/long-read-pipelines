version 1.0

import "../tasks/Utility/Utils.wdl" as Utils
import "../tasks/Utility/ONTUtils.wdl" as ONT
import "tasks/C3POa.wdl" as C3
import "../tasks/Alignment/AlignReads.wdl" as AR
import "../tasks/QC/AlignedMetrics.wdl" as AM
import "../tasks/Utility/Finalize.wdl" as FF

workflow ONT10x {
    input {
        Array[File] final_summaries
        Array[File] sequencing_summaries
        File ref_map_file

        String participant_name
        Int num_shards = 300

        String gcs_out_root_dir
    }

    parameter_meta {
        final_summaries:      "GCS path to '*final_summary*.txt*' files for basecalled fastq files"
        sequencing_summaries: "GCS path to '*sequencing_summary*.txt*' files for basecalled fastq files"
        ref_map_file:         "table indicating reference sequence and auxillary file locations"

        participant_name:     "name of the participant from whom these samples were obtained"
        num_shards:           "[default-valued] number of shards into which fastq files should be batched"

        gcs_out_root_dir:     "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONT10x/" + participant_name

    scatter (p in zip(final_summaries, sequencing_summaries)) {
        File final_summary = p.left
        File sequencing_summary = p.right

        call ONT.GetRunInfo { input: final_summary = final_summary }
        call ONT.ListFiles as ListFastqs { input: sequencing_summary = sequencing_summary, suffix = "fastq" }

        String SM  = participant_name
        String PL  = "ONT"
        String PU  = GetRunInfo.run_info["instrument"]
        String DT  = GetRunInfo.run_info["started"]
        String ID  = GetRunInfo.run_info["flow_cell_id"] + "." + GetRunInfo.run_info["position"]
        String DIR = GetRunInfo.run_info["protocol_group_id"] + "." + SM + "." + ID
        String SID = ID + "." + sub(GetRunInfo.run_info["protocol_run_id"], "-.*", "")
        String RG = "@RG\\tID:~{SID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        String rg_subreads  = "@RG\\tID:~{SID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = ListFastqs.manifest, N = num_shards }

        scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
#            call AR.Minimap2 as AlignSubreads {
#                input:
#                    reads      = read_lines(manifest_chunk),
#                    ref_fasta  = ref_map['fasta'],
#                    RG         = rg_subreads,
#                    map_preset = "splice"
#            }

            call C3.C3POa as C3POa { input: manifest_chunk = manifest_chunk, ref_fasta = ref_map['fasta'] }

            scatter (a in zip([1, 2, 3, 4], [ C3POa.consensus1, C3POa.consensus2, C3POa.consensus3, C3POa.consensus4 ])) {
                Int splint_num = a.left
                File fq = a.right

#                call FastaToSam { input: fasta = C3POa.consensus }
#                call AnnotateAdapters { input: bam = FastaToSam.output_bam }

                String rg_consensus = "@RG\\tID:~{SID}.consensus~{splint_num}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

                call AR.Minimap2 as AlignConsensus {
                    input:
                        reads      = [ fq ],
                        ref_fasta  = ref_map['fasta'],
                        RG         = rg_consensus,
                        map_preset = "splice"
                }

                #call CountFastaRecords as CountConsensusReadsInPartition { input: fasta = fq }
            }

#            File align_subreads_bam = AlignSubreads.aligned_bam

            File align_consensus_bam1 = AlignConsensus.aligned_bam[0]
            File align_consensus_bam2 = AlignConsensus.aligned_bam[1]
            File align_consensus_bam3 = AlignConsensus.aligned_bam[2]
            File align_consensus_bam4 = AlignConsensus.aligned_bam[3]

            call CountNumPasses as CountNumPasses1 { input: fastq = C3POa.subreads1 }
            call CountNumPasses as CountNumPasses2 { input: fastq = C3POa.subreads2 }
            call CountNumPasses as CountNumPasses3 { input: fastq = C3POa.subreads3 }
            call CountNumPasses as CountNumPasses4 { input: fastq = C3POa.subreads4 }

            call CountFastqRecords as CountSubreadsInPartition1 { input: fastq = C3POa.subreads1 }
            call CountFastqRecords as CountSubreadsInPartition2 { input: fastq = C3POa.subreads2 }
            call CountFastqRecords as CountSubreadsInPartition3 { input: fastq = C3POa.subreads3 }
            call CountFastqRecords as CountSubreadsInPartition4 { input: fastq = C3POa.subreads4 }

#            call CountFastqRecords as CountAnnotatedReadsInPartition { input: fastq = AnnotateAdapters.annotated_fq }

            call CountFastaRecords as CountConsensusReadsInPartition1 { input: fasta = C3POa.consensus1 }
            call CountFastaRecords as CountConsensusReadsInPartition2 { input: fasta = C3POa.consensus2 }
            call CountFastaRecords as CountConsensusReadsInPartition3 { input: fasta = C3POa.consensus3 }
            call CountFastaRecords as CountConsensusReadsInPartition4 { input: fasta = C3POa.consensus4 }
        }

        call Utils.Sum as CountNoSplintReadsInRun  { input: ints = C3POa.no_splint_reads }
        call Utils.Sum as CountUnderLenCutoffInRun { input: ints = C3POa.under_len_cutoff }
        call Utils.Sum as CountTotalSubreadsInRun  { input: ints = C3POa.total_reads }

        call C3.Cat as CountNumPassesInRun1 { input: files = CountNumPasses1.num_passes, out = "~{participant_name}.num_passes1.txt" }
        call C3.Cat as CountNumPassesInRun2 { input: files = CountNumPasses2.num_passes, out = "~{participant_name}.num_passes2.txt" }
        call C3.Cat as CountNumPassesInRun3 { input: files = CountNumPasses3.num_passes, out = "~{participant_name}.num_passes3.txt" }
        call C3.Cat as CountNumPassesInRun4 { input: files = CountNumPasses4.num_passes, out = "~{participant_name}.num_passes4.txt" }

        call Utils.Sum as CountSubreadsInRun1 { input: ints = CountSubreadsInPartition1.num_records }
        call Utils.Sum as CountSubreadsInRun2 { input: ints = CountSubreadsInPartition2.num_records }
        call Utils.Sum as CountSubreadsInRun3 { input: ints = CountSubreadsInPartition3.num_records }
        call Utils.Sum as CountSubreadsInRun4 { input: ints = CountSubreadsInPartition4.num_records }

#        call Utils.Sum as CountAnnotatedReadsInRun { input: ints = CountAnnotatedReadsInPartition.num_records }

        call Utils.Sum as CountConsensusReadsInRun1 { input: ints = CountConsensusReadsInPartition1.num_records }
        call Utils.Sum as CountConsensusReadsInRun2 { input: ints = CountConsensusReadsInPartition2.num_records }
        call Utils.Sum as CountConsensusReadsInRun3 { input: ints = CountConsensusReadsInPartition3.num_records }
        call Utils.Sum as CountConsensusReadsInRun4 { input: ints = CountConsensusReadsInPartition4.num_records }

#        call Utils.MergeBams as MergeAnnotated { input: bams = AnnotateAdapters.annotated_bam }
#        call Utils.MergeBams as MergeSubreads { input: bams = align_subreads_bam, prefix = "~{participant_name}.subreads" }

        call Utils.MergeBams as MergeConsensus1 { input: bams = align_consensus_bam1, prefix = "~{participant_name}.consensus1" }
        call Utils.MergeBams as MergeConsensus2 { input: bams = align_consensus_bam2, prefix = "~{participant_name}.consensus2" }
        call Utils.MergeBams as MergeConsensus3 { input: bams = align_consensus_bam3, prefix = "~{participant_name}.consensus3" }
        call Utils.MergeBams as MergeConsensus4 { input: bams = align_consensus_bam4, prefix = "~{participant_name}.consensus4" }
    }

    call Utils.Sum as CountNoSplintReads  { input: ints = CountNoSplintReadsInRun.sum,  prefix = "~{participant_name}.no_splint_reads" }
    call Utils.Sum as CountUnderLenCutoff { input: ints = CountUnderLenCutoffInRun.sum, prefix = "~{participant_name}.under_len_cutoff" }
    call Utils.Sum as CountTotalSubreads  { input: ints = CountTotalSubreadsInRun.sum,  prefix = "~{participant_name}.total_subreads" }

    call C3.Cat as CountNumPassesAll1 { input: files = CountNumPassesInRun1.merged, out = "~{participant_name}.num_passes1.txt" }
    call C3.Cat as CountNumPassesAll2 { input: files = CountNumPassesInRun2.merged, out = "~{participant_name}.num_passes2.txt" }
    call C3.Cat as CountNumPassesAll3 { input: files = CountNumPassesInRun3.merged, out = "~{participant_name}.num_passes3.txt" }
    call C3.Cat as CountNumPassesAll4 { input: files = CountNumPassesInRun4.merged, out = "~{participant_name}.num_passes4.txt" }

    call Utils.Sum as CountSubreads1 { input: ints = CountSubreadsInRun1.sum, prefix = "~{participant_name}.num_subreads1" }
    call Utils.Sum as CountSubreads2 { input: ints = CountSubreadsInRun2.sum, prefix = "~{participant_name}.num_subreads2" }
    call Utils.Sum as CountSubreads3 { input: ints = CountSubreadsInRun3.sum, prefix = "~{participant_name}.num_subreads3" }
    call Utils.Sum as CountSubreads4 { input: ints = CountSubreadsInRun4.sum, prefix = "~{participant_name}.num_subreads4" }

#    call Utils.Sum as CountAnnotatedReads { input: ints = CountAnnotatedReadsInRun.sum, prefix = "num_annotated" }

    call Utils.Sum as CountConsensusReads1 { input: ints = CountConsensusReadsInRun1.sum, prefix = "~{participant_name}.num_consensus1" }
    call Utils.Sum as CountConsensusReads2 { input: ints = CountConsensusReadsInRun2.sum, prefix = "~{participant_name}.num_consensus2" }
    call Utils.Sum as CountConsensusReads3 { input: ints = CountConsensusReadsInRun3.sum, prefix = "~{participant_name}.num_consensus3" }
    call Utils.Sum as CountConsensusReads4 { input: ints = CountConsensusReadsInRun4.sum, prefix = "~{participant_name}.num_consensus4" }

#    call Utils.MergeBams as MergeAllAnnotated { input: bams = MergeAnnotated.merged_bam, prefix = "~{participant_name}.annotated" }

    if (length(MergeConsensus1.merged_bam) > 1) {
#        call Utils.MergeBams as MergeAllSubreads { input: bams = MergeSubreads.merged_bam, prefix = "~{participant_name}.subreads" }

        call Utils.MergeBams as MergeAllConsensus1 { input: bams = MergeConsensus1.merged_bam, prefix = "~{participant_name}.consensus1" }
        call Utils.MergeBams as MergeAllConsensus2 { input: bams = MergeConsensus2.merged_bam, prefix = "~{participant_name}.consensus2" }
        call Utils.MergeBams as MergeAllConsensus3 { input: bams = MergeConsensus3.merged_bam, prefix = "~{participant_name}.consensus3" }
        call Utils.MergeBams as MergeAllConsensus4 { input: bams = MergeConsensus4.merged_bam, prefix = "~{participant_name}.consensus4" }
    }

#    File subreads_bam = select_first([MergeAllSubreads.merged_bam, MergeSubreads.merged_bam[0]])
#    File subreads_bai = select_first([MergeAllSubreads.merged_bai, MergeSubreads.merged_bai[0]])

    File consensus_bam1 = select_first([MergeAllConsensus1.merged_bam, MergeConsensus1.merged_bam[0]])
    File consensus_bai1 = select_first([MergeAllConsensus1.merged_bai, MergeConsensus1.merged_bai[0]])
    File consensus_bam2 = select_first([MergeAllConsensus2.merged_bam, MergeConsensus2.merged_bam[0]])
    File consensus_bai2 = select_first([MergeAllConsensus2.merged_bai, MergeConsensus2.merged_bai[0]])
    File consensus_bam3 = select_first([MergeAllConsensus3.merged_bam, MergeConsensus3.merged_bam[0]])
    File consensus_bai3 = select_first([MergeAllConsensus3.merged_bai, MergeConsensus3.merged_bai[0]])
    File consensus_bam4 = select_first([MergeAllConsensus4.merged_bam, MergeConsensus4.merged_bam[0]])
    File consensus_bai4 = select_first([MergeAllConsensus4.merged_bai, MergeConsensus4.merged_bai[0]])

#    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBC {
#        input:
#            bam = MergeAllAnnotated.merged_bam,
#            prefix = "num_annotated_with_cbc",
#            regex = "CB:Z:[ACGT]"
#    }
#
#    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBCAndUniqueAlignment {
#        input:
#            bam = MergeAllAnnotated.merged_bam,
#            samfilter = "-F 0x100",
#            prefix = "num_annotated_with_cbc_and_unique_alignment",
#            regex = "CB:Z:[ACGT]"
#    }
#
#    call Utils.BamToTable { input: bam = MergeAllAnnotated.merged_bam, prefix = "reads_aligned_annotated.table" }
#
#    call AM.AlignedMetrics as PerFlowcellRunConsensusMetrics {
#        input:
#            aligned_bam    = MergeAllConsensus.merged_bam,
#            aligned_bai    = MergeAllConsensus.merged_bai,
#            ref_fasta      = ref_map['fasta'],
#            ref_dict       = ref_map['dict'],
#            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
#    }
#
#    ##########
#    # Finalize
#    ##########

    call FF.FinalizeToDir as FinalizeReads {
        input:
            files = [ consensus_bam1, consensus_bai1,
                      consensus_bam2, consensus_bai2,
                      consensus_bam3, consensus_bai3,
                      consensus_bam4, consensus_bai4
                    ],
            outdir = outdir + "/alignments"
    }

    call FF.FinalizeToDir as FinalizeC3POaStats {
        input:
            files = [ CountNoSplintReads.sum_file, CountUnderLenCutoff.sum_file, CountTotalSubreads.sum_file ],
            outdir = outdir + "/metrics/c3poa_stats"
    }

    call FF.FinalizeToDir as FinalizeConsensusReadCounts {
        input:
            files = [ CountNumPassesAll1.merged, CountNumPassesAll2.merged, CountNumPassesAll3.merged, CountNumPassesAll4.merged,
                      CountSubreads1.sum_file, CountSubreads2.sum_file, CountSubreads3.sum_file, CountSubreads4.sum_file,
                      CountConsensusReads1.sum_file, CountConsensusReads2.sum_file, CountConsensusReads3.sum_file, CountConsensusReads4.sum_file
                    ],
            outdir = outdir + "/metrics/read_counts"
    }

#    call FF.FinalizeToDir as FinalizeConsensusReadCounts {
#        input:
#            files = [ CountSubreads.sum_file, CountAnnotatedReads.sum_file, CountConsensusReads.sum_file,
#                      GrepAnnotatedReadsWithCBC.num_records_file, GrepAnnotatedReadsWithCBCAndUniqueAlignment.num_records_file ],
#            outdir = outdir + "/metrics/read_counts"
#    }
#
#    call FF.FinalizeToDir as FinalizeNumPasses {
#        input:
#            files = [ CountNumPassesAll.merged ],
#            outdir = outdir + "/metrics/num_passes"
#    }
#
#    call FF.FinalizeToDir as FinalizeBamTable {
#        input:
#            files = [ BamToTable.table ],
#            outdir = outdir + "/metrics/bam_tables"
#    }
#
#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ MergeAllConsensus.merged_bam, MergeAllConsensus.merged_bai ],
#            outdir = outdir + "/alignments"
#    }
}

task AnnotateAdapters {
    input {
        File bam
        Int read_end_length = 500
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 2
    Int disk_size = 4*ceil(size(bam, "GB"))

    String output_name = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        python3 /lrma/tool.py \
            --bam=~{bam} \
            --adapter=/lrma/adapter_sequence.fasta \
            --reverse-adapter=/lrma/reverse_adapter_sequence.fasta \
            --whitelist-10x=/lrma/3M-february-2018.txt \
            --name=~{output_name}_annotated \
            --read-end-length=~{read_end_length} \
            --record-umis \
            --ssw-path /lrma/ssw/ \
            --starcode-path /lrma/starcode-master/starcode

        /opt/conda/envs/10x_tool/bin/samtools fastq -T ZA,CR,ZU,CB ~{output_name}_annotated.bam | gzip > ~{output_name}_annotated.fastq.gz
    >>>

    output {
        File annotated_bam  = "~{output_name}_annotated.bam"
        File annotated_fq   = "~{output_name}_annotated.fastq.gz"
        File barcode_stats  = "~{output_name}_annotated_barcode_stats.tsv"
        File starcode_stats = "~{output_name}_annotated_starcode.tsv"
        File stats          = "~{output_name}_annotated_stats.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.9"
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

task CountNumPasses {
    input {
        File fastq
        String prefix = "num_passes"
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 2
    Int disk_size = 2*ceil(size(fastq, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{fastq} | paste - - - - | awk '{ print $1 }' | sed 's/_[0-9]*$//' | uniq -c > ~{prefix}.txt
    >>>

    output {
        File num_passes = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
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

task FastaToSam {

    meta {
        description: "Convert a fasta file to a sam file"
    }

    parameter_meta {
        fasta: "The fasta file"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Float fasta_sam_disk_multiplier = 3.25
    Int disk_size = ceil(fasta_sam_disk_multiplier * size(fasta, "GiB")) + 20

    command <<<
        python /usr/local/bin/prepare_run.py ~{fasta}
    >>>

    output {
        File output_bam = "unmapped.bam"
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

task CountFastqRecords {

    meta {
        description: "Count the number of records in a fastq file"
    }

    parameter_meta {
        fastq: "The fastq file"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File fastq

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(fastq, "GiB"))

    command <<<
        set -euxo pipefail

        FILE="~{fastq}"
        if [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            cat ~{fastq} | awk '{s++}END{print s/4}'
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            zcat ~{fastq} | awk '{s++}END{print s/4}'
        fi
    >>>

    output {
        Int num_records = read_int(stdout())
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

task CountFastaRecords {

    meta {
        description: "Count the number of records in a fasta file"
    }

    parameter_meta {
        fasta: "The fasta file"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(fasta, "GiB"))

    command <<<
        grep -c '>' ~{fasta}

        exit 0
    >>>

    output {
        Int num_records = read_int(stdout())
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
