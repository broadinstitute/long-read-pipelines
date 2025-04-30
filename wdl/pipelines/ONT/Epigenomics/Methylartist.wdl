version 1.0

workflow methylartist{
    meta{
        description: "a workflow that using methylartist to plot methylated signal from MM/ML tagged bams"
    }

    input{
        Array[File] bams
        Array[File] bais
        File reference
        String region
        String motif
    }
    Int data_length = length(bams)
    Array[Int] indexes= range(data_length)
    scatter (idx in indexes) {
        call SubsetBam{
            input: bam = bams[idx],
                   bai = bais[idx],
                   locus = region,
                   prefix = "subset" + idx
        }
    }
    
    call Methylartist{input: bams=SubsetBam.subset_bam, bais=SubsetBam.subset_bai, reference = reference, region = region, motif = motif}
    
    output{
        Array[File] output_file = Methylartist.output_files
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


task Methylartist{
    input{
        Array[File] bams
        Array[File] bais
        File reference
        String region
        String motif = "CG"
    }

    command <<<
        methylartist locus -b ~{sep="," bams} -i ~{region} --ref ~{reference} --motif ~{motif}
    >>>

    output{
        Array[File] output_files = glob("*.png")
    }

    Int disk_size = 100 + ceil(2 * (size(bams, "GiB")))

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/methylartist:v1"
    }
}