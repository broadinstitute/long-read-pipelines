version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/HiFi.wdl" as HIFI
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/AlignedMetrics.wdl" as AM
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/UnalignedMetrics.wdl" as UM
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/AnnotateAdapters.wdl" as AA
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/Figures.wdl" as FIG
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.21/wdl/tasks/Finalize.wdl" as FF

workflow PB10xSingleProcessedSample {
    input {
        File input_bam
        String? sample_name

        Array[File] subread_bams

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

    scatter (subread_bam in subread_bams) {
        call GetRunInfo as GRI { input: input_bam = subread_bam }

        #call Utils.CountBamRecords as CountSubreadsInBamShard { input: bam = unmapped_shard }

        call UM.UnalignedMetrics as PerFlowcellUnalignedMetrics {
            input:
                unaligned_bam = subread_bam,

                per = "flowcell",
                type = "raw",
                label = GRI.run_info['ID'],

                gcs_output_dir = outdir + "/" + GRI.run_info['ID']
        }
    }

    call RevertBam { input: input_bam = input_bam }

    scatter (rg_bam in RevertBam.rg_bams) {
        call GetRunInfo { input: input_bam = rg_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String ID  = GetRunInfo.run_info["ID"]
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID

        String rg_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}"
        String rg_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}"

        call PB.ShardLongReads { input: unmapped_files = [ rg_bam ], num_reads_per_split = 100000 }

        scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
            call Utils.CountBamRecords as CountConsensusReadsInShard { input: bam = unmapped_shard }

            call AA.AnnotateAdapters { input: bam = unmapped_shard }

            call Utils.CountFastqRecords as CountAnnotatedReadsInShard { input: fastq = AnnotateAdapters.annotated_fq }

            call AR.Minimap2 as AlignSubreads {
                input:
                    reads = [ unmapped_shard ],
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

        call Utils.Sum as CountConsensusReads { input: ints = CountConsensusReadsInShard.num_records, prefix = "num_consensus" }
        call Utils.Sum as CountAnnotatedReads { input: ints = CountAnnotatedReadsInShard.num_records, prefix = "num_annotated" }

        call AR.MergeBams as MergeSubreads  { input: bams = AlignSubreads.aligned_bam, prefix = "~{SM}.~{ID}.subreads" }
        call AR.MergeBams as MergeConsensus { input: bams = AlignConsensus.aligned_bam, prefix = "~{SM}.~{ID}.consensus" }
        call AR.MergeBams as MergeAnnotated { input: bams = AnnotateAdapters.annotated_bam, prefix = "~{SM}.~{ID}.annotated" }

        call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBC {
            input:
                bam = MergeAnnotated.merged_bam,
                prefix = "num_annotated_with_cbc",
                regex = "CB:Z:[ACGT]"
        }

        call Utils.GrepCountBamRecords as GrepAnnotatedReadsWithCBCAndUniqueAlignment {
            input:
                bam = MergeAnnotated.merged_bam,
                samfilter = "-F 0x100",
                prefix = "num_annotated_with_cbc_and_unique_alignment",
                regex = "CB:Z:[ACGT]"
        }

        call Utils.BamToTable { input: bam = MergeAnnotated.merged_bam, prefix = "reads_aligned_annotated.table" }

        call AM.AlignedMetrics as PerFlowcellSubreadMetrics {
            input:
                aligned_bam    = MergeSubreads.merged_bam,
                aligned_bai    = MergeSubreads.merged_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                ref_flat       = ref_flat,
                dbsnp_vcf      = dbsnp_vcf,
                dbsnp_tbi      = dbsnp_tbi,
                metrics_locus  = metrics_locus,
                per            = "flowcell",
                type           = "subrun",
                label          = ID + ".subreads",
                gcs_output_dir = outdir + "/" + DIR
        }

        call AM.AlignedMetrics as PerFlowcellConsensusMetrics {
            input:
                aligned_bam    = MergeConsensus.merged_bam,
                aligned_bai    = MergeConsensus.merged_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                ref_flat       = ref_flat,
                dbsnp_vcf      = dbsnp_vcf,
                dbsnp_tbi      = dbsnp_tbi,
                metrics_locus  = metrics_locus,
                per            = "flowcell",
                type           = "subrun",
                label          = ID + ".consensus",
                gcs_output_dir = outdir + "/" + DIR
        }

        call AM.AlignedMetrics as PerFlowcellAnnotatedMetrics {
            input:
                aligned_bam    = MergeAnnotated.merged_bam,
                aligned_bai    = MergeAnnotated.merged_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                ref_flat       = ref_flat,
                dbsnp_vcf      = dbsnp_vcf,
                dbsnp_tbi      = dbsnp_tbi,
                metrics_locus  = metrics_locus,
                per            = "flowcell",
                type           = "run",
                label          = ID + ".annotated",
                gcs_output_dir = outdir + "/" + DIR
        }

        ##########
        # Finalize
        ##########

        call FF.FinalizeToDir as FinalizeReadCounts {
            input:
                files = [ CountConsensusReads.sum_file, CountAnnotatedReads.sum_file,
                          GrepAnnotatedReadsWithCBC.num_records_file,
                          GrepAnnotatedReadsWithCBCAndUniqueAlignment.num_records_file ],
                outdir = outdir + "/" + DIR + "/metrics/read_counts"
        }

        call FF.FinalizeToDir as FinalizeBamTable {
            input:
                files = [ BamToTable.table ],
                outdir = outdir + "/" + DIR + "/metrics/bam_tables"
        }

        call FF.FinalizeToDir as FinalizeMergedRuns {
            input:
                files = [ MergeConsensus.merged_bam, MergeConsensus.merged_bai,
                          MergeAnnotated.merged_bam, MergeAnnotated.merged_bai ],
                outdir = outdir + "/" + DIR + "/alignments"
        }

        call FF.FinalizeToDir as FinalizeMergedSubreads {
            input:
                files = [ MergeSubreads.merged_bam,  MergeSubreads.merged_bai ],
                outdir = outdir + "/" + DIR + "/alignments"
        }
    }
}

task RevertBam {
    input {
        File input_bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /usr/local/bin/gatk.jar RevertSam -I ~{input_bam} -OBR -O ./

        du -hcs .
    >>>

    output {
        Array[File] rg_bams = glob("*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             40,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
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

task SplitBam {
    input {
        File input_bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools split ~{input_bam}
    >>>

    output {
        Array[File] rg_bams = glob("*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
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

task GetRunInfo {
    input {
        String input_bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view -H ~{input_bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep -v '^@RG' | sed 's/:/\t/'
    >>>

    output {
        Map[String, String] run_info = read_map(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
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
