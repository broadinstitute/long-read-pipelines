version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Utility/cnv_common_tasks.wdl" as CNVTasks

import "MakeBincovMatrix.wdl" as MaMa
import "PloidyEstimation.wdl" as PE

workflow EstimateLrSamplePloidy {
  meta {
    description: "Use parts of the GATK-SV pipeline for sample ploidy estimation."
  }
  input {
    Array[File] bams
    Array[File] bais

    File reference_map
    File intervals

    ####################################################
    #### optional arguments for PreprocessIntervals ####
    ####################################################
    File? blacklist_intervals
    Int? padding

    File primary_contigs_list
    File primary_contigs_fai

    Int binsize
    String batch_name

    String gatk_docker
    String sv_pipeline_base_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_qc_docker
  }
  output {

  }

  Int mem_gb_for_collect_counts = 32
  Int disk_space_gb_for_collect_counts = 200

  Map[String, String] ref_map = read_map(reference_map)

  call CNVTasks.PreprocessIntervals {
    input:
      intervals = intervals,
      blacklist_intervals = blacklist_intervals,
      ref_fasta = ref_map['fasta'],
      ref_fasta_fai = ref_map['fai'],
      ref_fasta_dict = ref_map['dict'],
      padding = padding,
      bin_length = binsize,
      # gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      # preemptible_attempts = preemptible_attempts
    }
  # File preprocessed_intervals


  ######################################
  # collect counts for each sample
  scatter (idx in range(length(bams))) {
    File this_bam = bams[idx]
    File this_bai = bais[idx]

    call Utils.InferSampleName { input: bam = this_bam, bai = this_bai}
    String this_sample = InferSampleName.sample_name

    call CollectCounts {
      input:
        intervals = PreprocessIntervals.preprocessed_intervals,
        cram_or_bam = this_bam,
        cram_or_bam_idx = this_bai,
        sample_id = this_sample,
        ref_fasta = ref_map['fasta'],
        ref_fasta_fai = ref_map['fai'],
        ref_fasta_dict = ref_map['dict'],
        gatk_docker = gatk_docker,
        mem_gb = mem_gb_for_collect_counts,
        disk_space_gb = disk_space_gb_for_collect_counts,
        disabled_read_filters = ["MappingQualityReadFilter"]
    }

    File this_read_counts = CollectCounts.counts

    call CountsMetrics {
      input:
        counts_file = this_read_counts,
        sample_id = this_sample,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  Array[String] all_sample_names = this_sample
  Array[File] all_read_counts = this_read_counts

  ######################################
  # make matrix then QC
  call MaMa.MakeBincovMatrix {
    input:
      samples = all_sample_names,
      count_files = all_read_counts,
      batch = batch_name,
      binsize = binsize,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker
  }
  call PE.Ploidy {
    input:
      bincov_matrix = MakeBincovMatrix.merged_bincov,
      batch = batch_name,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker
  }
}

task CollectCounts {
  input {
    File intervals
    File cram_or_bam
    File cram_or_bam_idx
    String sample_id
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File? gatk4_jar_override
    Array[String]? disabled_read_filters

    # Runtime parameters
    String gatk_docker
    Float? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts
  }

  parameter_meta {
    cram_or_bam: {
      localization_optional: true
    }
    cram_or_bam_idx: {
      localization_optional: true
    }
  }

  Float mem_overhead_gb = 2.0
  Float machine_mem_gb = select_first([mem_gb, 12.0])
  Int command_mem_mb = floor((machine_mem_gb - mem_overhead_gb) * 1024)
  Array[String] disabled_read_filters_arr = if(defined(disabled_read_filters))
    then
      prefix(
        "--disable-read-filter ",
        select_first([disabled_read_filters])
      )
    else
      []

  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

    gatk --java-options "-Xmx~{command_mem_mb}m" CollectReadCounts \
      -L ~{intervals} \
      --input ~{cram_or_bam} \
      --read-index ~{cram_or_bam_idx} \
      --reference ~{ref_fasta} \
      --format TSV \
      --interval-merging-rule OVERLAPPING_ONLY \
      --output ~{sample_id}.counts.tsv \
      ~{sep=' ' disabled_read_filters_arr}

    sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:~{sample_id}/g" ~{sample_id}.counts.tsv
    bgzip ~{sample_id}.counts.tsv
  >>>

  runtime {
    docker: gatk_docker
    memory: machine_mem_gb + " GiB"
    disks: "local-disk " + select_first([disk_space_gb, 10]) + if use_ssd then " SSD" else " HDD"
    cpu: select_first([cpu, 1])
    preemptible: select_first([preemptible_attempts, 3])
    maxRetries: 1
  }

  output {
    File counts = "~{sample_id}.counts.tsv.gz"
  }
}

task CountsMetrics {
  input {
    File counts_file
    String sample_id
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{sample_id}.raw-counts.tsv"
  }
  command <<<

    svtest raw-counts ~{counts_file} ~{sample_id} > ~{sample_id}.raw-counts.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_base_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
}
