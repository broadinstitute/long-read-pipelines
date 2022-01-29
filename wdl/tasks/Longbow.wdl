version 1.0

import "Structs.wdl"
import "Utils.wdl"

workflow Process {
    input {
        File bam

        String prefix = "out"
        String? model
        File? barcode_tag
        File? barcode_allowlist
        Boolean same_barcode_per_read = false
    }

    if (!defined(model)) {
        call Peek { input: bam = bam, n = 1000 }
    }

    String lbmodel = select_first([Peek.model, "mas15v2"])

    call Annotate { input: bam = bam, model = lbmodel }

    if (defined(barcode_tag) && defined(barcode_allowlist)) {
        call Refine {
            input:
                bam = Annotate.annotated_bam,
                barcode_tag = select_first([barcode_tag]),
                barcode_allowlist = select_first([barcode_allowlist]),
                same_barcode_per_read = same_barcode_per_read
        }
    }

    call Extract { input: bam = select_first([Refine.refined_bam, Annotate.annotated_bam]) }

    output {
        File extracted_bam = Extract.extracted_bam
    }
}

task Peek {
    input {
        File bam
        Int n = 1000

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: {
                 localization_optional: true
             }
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        ((samtools view -H ~{bam}) && (samtools view ~{bam} | head -n ~{n})) | samtools view -b > subset.bam

        longbow peek -n ~{n} -o model.txt subset.bam
    >>>

    output {
        String model = read_string("model.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:kvg_docker_pathfix"
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

task Annotate {
    input {
        File bam
        String model

        String prefix = "out"
        Int num_cpus = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow annotate -m ~{model} -t ~{num_cpus} -o ~{prefix}.annotated.bam ~{bam}
    >>>

    output {
        File annotated_bam = "~{prefix}.annotated.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:kvg_docker_pathfix"
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

task Refine {
    input {
        File bam

        File barcode_tag
        File barcode_allowlist
        Boolean same_barcode_per_read

        String prefix = "out"
        Int num_cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow refine \
            -t ~{num_cpus} \
            -b ~{barcode_tag} \
            -a ~{barcode_allowlist} \
            ~{true='-s' false='' same_barcode_per_read} \
            -o ~{prefix}.refined.bam \
            ~{bam}
    >>>

    output {
        File refined_bam = "~{prefix}.refined.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:kvg_docker_pathfix"
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

task Extract {
    input {
        File bam

        Int num_cpus = 2
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow filter ~{bam} | longbow segment | longbow extract -o ~{prefix}.extracted.bam
    >>>

    output {
        File extracted_bam = "~{prefix}.extracted.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:kvg_docker_pathfix"
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

task Stats {
    input {
        File bam

        String prefix = "stats"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow stats -o ~{prefix} ~{bam}
    >>>

    output {
        Array[File] pngs = glob("*.png")
        Array[File] svgs = glob("*.svg")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:kvg_docker_pathfix"
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
