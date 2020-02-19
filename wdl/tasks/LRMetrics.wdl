version 1.0

import "Utils.wdl" as Utils
import "AlignedMetrics.wdl" as AM
import "UnalignedMetrics.wdl" as UM
import "Figures.wdl" as FIG
import "Finalize.wdl" as FF

workflow LRMetrics {
    input {
        File bam
        File bai

        String sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        File ref_flat

        File dbsnp
        File dbsnp_tbi
        File metrics_locus

        String gcs_output_dir
    }

    String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

    call SplitByReadGroup { input: bam = bam }

    scatter (bam in SplitByReadGroup.bams) {
        call IndexBam { input: input_bam = bam }

#        call FIG.Figures {
#            input:
#                summary_file = PrepareData.summary_file,
#
#                per = "flowcell",
#                type = "raw",
#                label = PrepareData.run_info['ID'],
#
#                gcs_output_dir = gcs_output_dir
#        }

#        call UM.UnalignedMetrics as PerFlowcellUnalignedConsensusMetrics {
#            input:
#                unaligned_bam = IndexBam.bam,
#
#                per = "flowcell",
#                type = "raw",
#                label = basename(IndexBam.bam, ".bam"),
#
#                gcs_output_dir = gcs_output_dir
#        }

        call AM.AlignedMetrics as PerFlowcellConsensusMetrics {
            input:
                aligned_bam = IndexBam.bam,
                aligned_bai = IndexBam.bai,
                ref_fasta = ref_fasta,
                ref_dict = ref_dict,
                ref_flat = ref_flat,
                dbsnp_vcf = dbsnp,
                dbsnp_tbi = dbsnp_tbi,
                metrics_locus = metrics_locus,
                per = "flowcell",
                type = "consensus",
                label = basename(IndexBam.bam, ".bam"),
                gcs_output_dir = gcs_output_dir
        }

#        call AM.AlignedMetrics as PerFlowcellSubreadsMetrics {
#            input:
#                aligned_bam = C3POa.subreads_bam,
#                aligned_bai = C3POa.subreads_bai,
#                ref_fasta = ref_fasta,
#                ref_dict = ref_dict,
#                ref_flat = ref_flat,
#                dbsnp = dbsnp,
#                dbsnp_tbi = dbsnp_tbi,
#                metrics_locus = metrics_locus,
#                per = "flowcell",
#                type = "subreads",
#                label = PrepareData.run_info['ID'],
#                gcs_output_dir = gcs_output_dir
#        }
    }

    call AM.AlignedMetrics as PerSampleConsensusMetrics {
        input:
            aligned_bam = bam,
            aligned_bai = bai,
            ref_fasta = ref_fasta,
            ref_dict = ref_dict,
            ref_flat = ref_flat,
            dbsnp_vcf = dbsnp,
            dbsnp_tbi = dbsnp_tbi,
            metrics_locus = metrics_locus,
            per = "sample",
            type = "consensus",
            label = sample_name,
            gcs_output_dir = gcs_output_dir
    }

#    call AM.AlignedMetrics as PerSampleSubreadsMetrics {
#        input:
#            aligned_bam = MergeFinalSubreads.merged,
#            aligned_bai = MergeFinalSubreads.merged_bai,
#            ref_fasta = ref_fasta,
#            ref_dict = ref_dict,
#            ref_flat = ref_flat,
#            dbsnp = dbsnp,
#            dbsnp_tbi = dbsnp_tbi,
#            metrics_locus = metrics_locus,
#            per = "sample",
#            type = "subreads",
#            label = sample_name,
#            gcs_output_dir = gcs_output_dir
#    }
}

task SplitByReadGroup {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools split -f '%!.%.' ~{bam}
    >>>

    output {
        Array[File] bams = glob("*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.22"
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

task IndexBam {
    input {
        File input_bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(input_bam)

    Int disk_size = 2*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        mv ~{input_bam} ~{basename}
        samtools index ~{basename}
    >>>

    output {
        File bam = "~{basename}"
        File bai = "~{basename}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.22"
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
