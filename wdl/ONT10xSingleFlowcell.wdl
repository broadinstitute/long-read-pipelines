version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/ONTUtils.wdl" as ONT
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/C3POa.wdl" as C3
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/AlignedMetrics.wdl" as AM
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Figures.wdl" as FIG
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.27/wdl/tasks/Finalize.wdl" as FF

workflow ONT10xSingleFlowcell {
    input {
        String gcs_input_dir
        String? sample_name
        Int fastq_shards = 50

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File ref_flat
        File dbsnp_vcf
        File dbsnp_tbi

        File metrics_locus

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call ONT.FindSequencingSummaryFiles { input: gcs_input_dir = gcs_input_dir }

    scatter (summary_file in FindSequencingSummaryFiles.summary_files) {
        call ONT.GetRunInfo { input: summary_file = summary_file }
        #call ONT.ListFiles as ListFast5s { input: summary_file = summary_file, suffix = "fast5" }
        call ONT.ListFiles as ListFastqs { input: summary_file = summary_file, suffix = "fastq" }

        String SM  = select_first([sample_name, GetRunInfo.run_info["sample_id"]])
        String PL  = "ONT"
        String PU  = GetRunInfo.run_info["instrument"]
        String DT  = GetRunInfo.run_info["started"]
        String ID  = GetRunInfo.run_info["flow_cell_id"] + "." + GetRunInfo.run_info["position"]
        String DIR = GetRunInfo.run_info["protocol_group_id"] + "." + SM + "." + ID
        String SID = ID + "." + sub(GetRunInfo.run_info["protocol_run_id"], "-.*", "")

        String rg_subreads  = "@RG\\tID:~{SID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String rg_consensus = "@RG\\tID:~{SID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        #call ONT.PartitionManifest as PartitionFast5Manifest { input: manifest = ListFast5s.manifest, N = 4  }
        call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = ListFastqs.manifest, N = fastq_shards }

        scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
            call C3.C3POa as C3POa { input: manifest_chunk = manifest_chunk, ref_fasta = ref_fasta }

            call Utils.CountFastqRecords as CountSubreadsInPartition { input: fastq = C3POa.subreads }

            call Utils.FastaToSam as FastaToSam { input: fasta = C3POa.consensus }

            call Utils.CountFastaRecords as CountConsensusReadsInPartition { input: fasta = C3POa.consensus }

            call AnnotateAdapters { input: bam = FastaToSam.output_bam }

            call Utils.CountFastqRecords as CountAnnotatedReadsInPartition { input: fastq = AnnotateAdapters.annotated_fq }

            call AR.Minimap2 as AlignSubreads {
                input:
                    reads = [ C3POa.subreads ],
                    ref_fasta = ref_fasta,
                    RG = rg_subreads,
                    map_preset = "splice"
            }

            call AR.Minimap2 as AlignConsensus {
                input:
                    reads = [ AnnotateAdapters.annotated_fq ],
                    ref_fasta = ref_fasta,
                    RG = rg_consensus,
                    map_preset = "splice"
            }
        }

        call Utils.Sum as CountSubreadsInRun { input: ints = CountSubreadsInPartition.num_records }
        call Utils.Sum as CountConsensusReadsInRun { input: ints = CountConsensusReadsInPartition.num_records }
        call Utils.Sum as CountAnnotatedReadsInRun { input: ints = CountAnnotatedReadsInPartition.num_records }

        call Utils.MergeBams as MergeSubreads  { input: bams = AlignSubreads.aligned_bam }
        call Utils.MergeBams as MergeConsensus { input: bams = AlignConsensus.aligned_bam }
        call Utils.MergeBams as MergeAnnotated { input: bams = AnnotateAdapters.annotated_bam }

#        call AM.AlignedMetrics as PerFlowcellSubRunSubreadMetrics {
#            input:
#                aligned_bam    = MergeSubreads.merged_bam,
#                aligned_bai    = MergeSubreads.merged_bai,
#                ref_fasta      = ref_fasta,
#                ref_dict       = ref_dict,
#                ref_flat       = ref_flat,
#                dbsnp_vcf      = dbsnp_vcf,
#                dbsnp_tbi      = dbsnp_tbi,
#                metrics_locus  = metrics_locus,
#                per            = "flowcell",
#                type           = "subrun",
#                label          = SID + ".subreads",
#                gcs_output_dir = outdir + "/" + DIR
#        }

#        call AM.AlignedMetrics as PerFlowcellSubRunConsensusMetrics {
#            input:
#                aligned_bam    = MergeConsensus.merged_bam,
#                aligned_bai    = MergeConsensus.merged_bai,
#                ref_fasta      = ref_fasta,
#                ref_dict       = ref_dict,
#                ref_flat       = ref_flat,
#                dbsnp_vcf      = dbsnp_vcf,
#                dbsnp_tbi      = dbsnp_tbi,
#                metrics_locus  = metrics_locus,
#                per            = "flowcell",
#                type           = "subrun",
#                label          = SID + ".consensus",
#                gcs_output_dir = outdir + "/" + DIR
#        }

        call FIG.Figures as PerFlowcellSubRunFigures {
            input:
                summary_files  = [ summary_file ],

                per            = "flowcell",
                type           = "subrun",
                label          = SID,

                gcs_output_dir = outdir + "/" + DIR
        }
    }

    call Utils.Sum as CountSubreads { input: ints = CountSubreadsInRun.sum, prefix = "num_subreads" }
    call Utils.Sum as CountConsensusReads { input: ints = CountConsensusReadsInRun.sum, prefix = "num_consensus" }
    call Utils.Sum as CountAnnotatedReads { input: ints = CountAnnotatedReadsInRun.sum, prefix = "num_annotated" }

    #call Utils.MergeBams as MergeAllSubreads  { input: bams = MergeSubreads.merged_bam,  prefix = "~{SM[0]}.~{ID[0]}.subreads"  }
    call Utils.MergeBams as MergeAllConsensus { input: bams = MergeConsensus.merged_bam, prefix = "~{SM[0]}.~{ID[0]}.consensus" }
    call Utils.MergeBams as MergeAllAnnotated { input: bams = MergeAnnotated.merged_bam, prefix = "~{SM[0]}.~{ID[0]}.annotated" }

#    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBC {
#        input:
#            bam = MergeAllAnnotated.merged_bam,
#            prefix = "num_annotated_with_cbc",
#            regex = "CB:Z:[ACGT]"
#    }

#    call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBCAndUniqueAlignment {
#        input:
#            bam = MergeAllAnnotated.merged_bam,
#            samfilter = "-F 0x100",
#            prefix = "num_annotated_with_cbc_and_unique_alignment",
#            regex = "CB:Z:[ACGT]"
#    }

    call Utils.BamToTable { input: bam = MergeAllAnnotated.merged_bam, prefix = "reads_aligned_annotated.table" }

#    call AM.AlignedMetrics as PerFlowcellRunSubreadMetrics {
#        input:
#            aligned_bam    = MergeAllSubreads.merged_bam,
#            aligned_bai    = MergeAllSubreads.merged_bai,
#            ref_fasta      = ref_fasta,
#            ref_dict       = ref_dict,
#            ref_flat       = ref_flat,
#            dbsnp_vcf      = dbsnp_vcf,
#            dbsnp_tbi      = dbsnp_tbi,
#            metrics_locus  = metrics_locus,
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0] + ".subreads",
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }

#    call AM.AlignedMetrics as PerFlowcellRunConsensusMetrics {
#        input:
#            aligned_bam    = MergeAllConsensus.merged_bam,
#            aligned_bai    = MergeAllConsensus.merged_bai,
#            ref_fasta      = ref_fasta,
#            ref_dict       = ref_dict,
#            ref_flat       = ref_flat,
#            dbsnp_vcf      = dbsnp_vcf,
#            dbsnp_tbi      = dbsnp_tbi,
#            metrics_locus  = metrics_locus,
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0] + ".consensus",
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }

#    call AM.AlignedMetrics as PerFlowcellRunAnnotatedMetrics {
#        input:
#            aligned_bam    = MergeAllAnnotated.merged_bam,
#            aligned_bai    = MergeAllAnnotated.merged_bai,
#            ref_fasta      = ref_fasta,
#            ref_dict       = ref_dict,
#            ref_flat       = ref_flat,
#            dbsnp_vcf      = dbsnp_vcf,
#            dbsnp_tbi      = dbsnp_tbi,
#            metrics_locus  = metrics_locus,
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0] + ".annotated",
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }

    call FIG.Figures as PerFlowcellRunFigures {
        input:
            summary_files  = FindSequencingSummaryFiles.summary_files,

            per            = "flowcell",
            type           = "run",
            label          = ID[0],

            gcs_output_dir = outdir + "/" + DIR[0]
    }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeReadCounts {
        input:
            files = [ CountSubreads.sum_file, CountConsensusReads.sum_file, CountAnnotatedReads.sum_file ],
                      #GrepAnnotatedReadsWithCBC.num_records_file,
                      #GrepAnnotatedReadsWithCBCAndUniqueAlignment.num_records_file ],
            outdir = outdir + "/" + DIR[0] + "/metrics/read_counts"
    }

    call FF.FinalizeToDir as FinalizeBamTable {
        input:
            files = [ BamToTable.table ],
            outdir = outdir + "/" + DIR[0] + "/metrics/bam_tables"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ MergeAllConsensus.merged_bam, MergeAllConsensus.merged_bai,
                      MergeAllAnnotated.merged_bam, MergeAllAnnotated.merged_bai ],
            outdir = outdir + "/" + DIR[0] + "/alignments"
    }

#    call FF.FinalizeToDir as FinalizeMergedSubreads {
#        input:
#            files = [ MergeAllSubreads.merged_bam,  MergeAllSubreads.merged_bai ],
#            outdir = outdir + "/" + DIR[0] + "/alignments"
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
