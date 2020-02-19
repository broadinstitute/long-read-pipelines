version 1.0

import "Structs.wdl"
import "ProcessReads.wdl" as PR
import "Finalize.wdl" as FF

workflow SRWholeGenomeSingleSample {
    input {
        Array[String] gcs_dirs

        String sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_pac
        File ref_fasta_sa

        String gcs_output_dir
    }

    String outdir = sub(sub(gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")

    scatter (gcs_dir in gcs_dirs) {
        call GetFileList { input: gcs_dir = gcs_dir }

        call BwaMem as BwaMemPairedEnd {
            input:
                fastqs        = GetFileList.pairs,

                SM            = sample_name,
                ID            = basename(gcs_dir),

                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_dict      = ref_dict,
                ref_fasta_amb = ref_fasta_amb,
                ref_fasta_ann = ref_fasta_ann,
                ref_fasta_bwt = ref_fasta_bwt,
                ref_fasta_pac = ref_fasta_pac,
                ref_fasta_sa  = ref_fasta_sa
        }

        call BwaMem as BwaMemSingleEnd {
            input:
                fastqs        = [ GetFileList.single ],

                SM            = sample_name,
                ID            = basename(gcs_dir),

                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_dict      = ref_dict,
                ref_fasta_amb = ref_fasta_amb,
                ref_fasta_ann = ref_fasta_ann,
                ref_fasta_bwt = ref_fasta_bwt,
                ref_fasta_pac = ref_fasta_pac,
                ref_fasta_sa  = ref_fasta_sa
        }

        call PR.MergeBams as MergeShards {
            input:
                shards = [ BwaMemPairedEnd.bam, BwaMemSingleEnd.bam ],
                merged_name = "merged_shards.bam"
        }
    }

    call PR.MergeBams as MergeRuns {
        input:
            shards = MergeShards.merged,
            merged_name = "~{sample_name}.bam"
    }

    call FF.FinalizeToDir as FinalizeBams {
        input:
            files = [ MergeRuns.merged, MergeRuns.merged_bai ],
            outdir = outdir + "/alignments"
    }
}

task GetFileList {
    input {
        String gcs_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(sub(gcs_dir + "/", "/+", "/"), "gs:/", "gs://")

    command <<<
        set -euxo pipefail

        gsutil ls ~{indir}**_[12].fastq.gz | head -2 > pairs.txt
        gsutil ls ~{indir}**.fastq.gz | grep -v '_' | head -1 > singles.txt
    >>>

    output {
        Array[String] pairs = read_lines("pairs.txt")
        String single = read_string("singles.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-utils:0.01.03"
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

task BwaMem {
    input {
        Array[File] fastqs

        String SM
        String ID

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_pac
        File ref_fasta_sa

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 16
    Int disk_size = 6*ceil(size(fastqs, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB") + size(ref_fasta_amb, "GB") + size(ref_fasta_ann, "GB") + size(ref_fasta_bwt, "GB") + size(ref_fasta_pac, "GB") + size(ref_fasta_sa, "GB"))

    command <<<
        set -euxo pipefail

        bwa mem -t ~{cpus} -R "@RG\tID:~{ID}\tSM:~{SM}\tPL:ILLUMINA" ~{ref_fasta} ~{sep=' ' fastqs} > aligned.sam
        samtools sort -@~{cpus} -m4G -o aligned.bam aligned.sam
        samtools index -@2 aligned.bam
    >>>

    output {
        File sam = "aligned.sam"
        File bam = "aligned.bam"
        File bai = "aligned.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             60,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/sr-utils:0.01.00"
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
