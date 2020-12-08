version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/ONTUtils.wdl" as ONT
import "tasks/C3POa.wdl" as C3
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF

workflow ONT10x {
    input {
        Array[File] final_summaries
        Array[File] sequencing_summaries
        File ref_map_file

        String participant_name
        Int num_shards = 50

        String? gcs_out_root_dir
    }

    parameter_meta {
        final_summaries:           "GCS path to '*final_summary*.txt*' files for basecalled fastq files"
        sequencing_summaries:      "GCS path to '*sequencing_summary*.txt*' files for basecalled fastq files"
        ref_map_file:              "table indicating reference sequence and auxillary file locations"

        participant_name:          "name of the participant from whom these samples were obtained"
        num_shards:                "[default-valued] number of shards into which fastq files should be batched"

        gcs_out_root_dir:          "[optional] GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call Utils.GetDefaultDir { input: workflow_name = "ONT10x" }
    String outdir = sub(select_first([gcs_out_root_dir, GetDefaultDir.path]), "/$", "") + "/" + participant_name

    scatter (p in zip(final_summaries, sequencing_summaries)) {
        File final_summary = p.left
        File sequencing_summary = p.right

        call ONT.GetRunInfo { input: summary_file = final_summary }
        call ONT.ListFiles as ListFastqs { input: summary_file = final_summary, suffix = "fastq" }

        String SM  = participant_name
        String PL  = "ONT"
        String PU  = GetRunInfo.run_info["instrument"]
        String DT  = GetRunInfo.run_info["started"]
        String ID  = GetRunInfo.run_info["flow_cell_id"] + "." + GetRunInfo.run_info["position"]
        String DIR = GetRunInfo.run_info["protocol_group_id"] + "." + SM + "." + ID
        String SID = ID + "." + sub(GetRunInfo.run_info["protocol_run_id"], "-.*", "")
        String RG = "@RG\\tID:~{SID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        String rg_subreads  = "@RG\\tID:~{SID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String rg_consensus = "@RG\\tID:~{SID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = ListFastqs.manifest, N = num_shards }

        scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
            call C3.C3POa as C3POa { input: manifest_chunk = manifest_chunk, ref_fasta = ref_map['fasta'] }

            call Utils.FastaToSam as FastaToSam { input: fasta = C3POa.consensus }
            call AnnotateAdapters { input: bam = FastaToSam.output_bam }

            call AR.Minimap2 as AlignConsensus {
                input:
                    reads      = [ AnnotateAdapters.annotated_fq ],
                    ref_fasta  = ref_map['fasta'],
                    RG         = rg_consensus,
                    map_preset = "splice"
            }

            call CountNumPasses { input: fastq = C3POa.subreads }

            call Utils.CountFastqRecords as CountSubreadsInPartition { input: fastq = C3POa.subreads }
            call Utils.CountFastqRecords as CountAnnotatedReadsInPartition { input: fastq = AnnotateAdapters.annotated_fq }
            call Utils.CountFastaRecords as CountConsensusReadsInPartition { input: fasta = C3POa.consensus }
        }

        call C3.Cat as CountNumPassesInRun { input: files = CountNumPasses.num_passes, out = "num_passes.txt" }

        call Utils.Sum as CountSubreadsInRun { input: ints = CountSubreadsInPartition.num_records }
        call Utils.Sum as CountAnnotatedReadsInRun { input: ints = CountAnnotatedReadsInPartition.num_records }
        call Utils.Sum as CountConsensusReadsInRun { input: ints = CountConsensusReadsInPartition.num_records }

        call Utils.MergeBams as MergeAnnotated { input: bams = AnnotateAdapters.annotated_bam }
        call Utils.MergeBams as MergeConsensus { input: bams = AlignConsensus.aligned_bam }

        call FIG.Figures as PerFlowcellSubRunFigures {
            input:
                summary_files  = [ sequencing_summary ],
                gcs_output_dir = outdir + "/" + DIR
        }
    }

    call C3.Cat as CountNumPassesAll { input: files = CountNumPassesInRun.merged, out = "num_passes.txt" }

    call Utils.Sum as CountSubreads { input: ints = CountSubreadsInRun.sum, prefix = "num_subreads" }
    call Utils.Sum as CountAnnotatedReads { input: ints = CountAnnotatedReadsInRun.sum, prefix = "num_annotated" }
    call Utils.Sum as CountConsensusReads { input: ints = CountConsensusReadsInRun.sum, prefix = "num_consensus" }

    call Utils.MergeBams as MergeAllAnnotated { input: bams = MergeAnnotated.merged_bam, prefix = "~{participant_name}.annotated" }
    call Utils.MergeBams as MergeAllConsensus { input: bams = MergeConsensus.merged_bam, prefix = "~{participant_name}.consensus" }

    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBC {
        input:
            bam = MergeAllAnnotated.merged_bam,
            prefix = "num_annotated_with_cbc",
            regex = "CB:Z:[ACGT]"
    }

    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBCAndUniqueAlignment {
        input:
            bam = MergeAllAnnotated.merged_bam,
            samfilter = "-F 0x100",
            prefix = "num_annotated_with_cbc_and_unique_alignment",
            regex = "CB:Z:[ACGT]"
    }

    call Utils.BamToTable { input: bam = MergeAllAnnotated.merged_bam, prefix = "reads_aligned_annotated.table" }

    call AM.AlignedMetrics as PerFlowcellRunConsensusMetrics {
        input:
            aligned_bam    = MergeAllConsensus.merged_bam,
            aligned_bai    = MergeAllConsensus.merged_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            ref_flat       = ref_map['flat'],
            dbsnp_vcf      = ref_map['dbsnp_vcf'],
            dbsnp_tbi      = ref_map['dbsnp_tbi'],
            metrics_locus  = ref_map['metrics_locus'],
            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

    call FIG.Figures as PerFlowcellRunFigures {
        input:
            summary_files  = sequencing_summaries,
            gcs_output_dir = outdir + "/" + DIR[0]
    }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeConsensusReadCounts {
        input:
            files = [ CountSubreads.sum_file, CountAnnotatedReads.sum_file, CountConsensusReads.sum_file,
                      GrepAnnotatedReadsWithCBC.num_records_file, GrepAnnotatedReadsWithCBCAndUniqueAlignment.num_records_file ],
            outdir = outdir + "/metrics/read_counts"
    }

    call FF.FinalizeToDir as FinalizeNumPasses {
        input:
            files = [ CountNumPassesAll.merged ],
            outdir = outdir + "/metrics/num_passes"
    }

    call FF.FinalizeToDir as FinalizeBamTable {
        input:
            files = [ BamToTable.table ],
            outdir = outdir + "/metrics/bam_tables"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ MergeAllConsensus.merged_bam, MergeAllConsensus.merged_bai ],
            outdir = outdir + "/alignments"
    }
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
        boot_disk_gb:       10,
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
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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
