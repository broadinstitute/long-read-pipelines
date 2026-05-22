version 1.0

workflow pbCpG{
    meta{
        description: "a workflow that using pbCpG to extract methylated signal from MM/ML tagged bams"
    }

    input{
        File bam
        File bai
        File reference_fa
        String locus
        String output_prefix
        Int nthreads
    }

    call PbCpG{input: bam=bam, bai=bai, reference_fa = reference_fa, outputprefix = output_prefix, num_threads = nthreads}

    output{
        File combined_bed = PbCpG.combined_bed
        File combined_bed_index = PbCpG.combined_bed_index
        File combined_bed_bw = PbCpG.combined_bed_bw
        File hap1_bed = PbCpG.hap1_bed
        File hap1_bed_index = PbCpG.hap1_bed_index
        File hap1_bed_bw = PbCpG.hap1_bed_bw
        File hap2_bed = PbCpG.hap2_bed
        File hap2_bed_index = PbCpG.hap2_bed_index
        File hap2_bed_bw = PbCpG.hap2_bed_bw
        File log = PbCpG.log
    }
}

task PbCpG{
    input{
        File bam
        File bai
        File reference_fa
        String mode = "model"
        String outputprefix
        Int num_threads
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        /pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
            --bam ~{bam} \
            --pileup-mode ~{mode} \
            --modsites-mode "denovo" \
            --ref ~{reference_fa} \
            --output-prefix ~{outputprefix} \
            --threads ~{num_threads}
    >>>

    output{
        File combined_bed = "~{outputprefix}.combined.bed.gz"
        File combined_bed_index = "~{outputprefix}.combined.bed.gz.tbi"
        File combined_bed_bw = "~{outputprefix}.combined.bw"
        File hap1_bed = "~{outputprefix}.hap1.bed.gz"
        File hap1_bed_index = "~{outputprefix}.hap1.bed.gz.tbi"
        File hap1_bed_bw = "~{outputprefix}.hap1.bw"
        File hap2_bed = "~{outputprefix}.hap2.bed.gz"
        File hap2_bed_index = "~{outputprefix}.hap2.bed.gz.tbi"
        File hap2_bed_bw = "~{outputprefix}.hap2.bw"
        File log = "~{outputprefix}.log"
    }

    Int disk_size = 100 + ceil(2 * (size(bam, "GiB")))

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/pbcpg:v3.0.0"
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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