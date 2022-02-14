version 1.0

import "Structs.wdl"
import "Utils.wdl"

workflow Process {
    input {
        File bam

        String prefix = "out"
        String model = "mas15v2"
        File? barcode_tag
        File? barcode_allowlist
        Boolean same_barcode_per_read = false
    }

    if (!defined(model)) { call Peek { input: bam = bam } }
    String lbmodel = select_first([model, Peek.model])

    call Annotate { input: bam = bam, model = lbmodel }
    call Filter   { input: bam = Annotate.annotated_bam }
    call Extract  { input: bam = Filter.filtered_bam }
    call Correct  { input: bam = Filter.filtered_bam, model = lbmodel }

    output {
        File annotated_bam = Annotate.annotated_bam
        File filtered_bam = Filter.filtered_bam
        File corrected_bam = Correct.corrected_bam
    }
}

task Peek {
    input {
        File bam
        Int n = 100

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow peek -n ~{n} -o model.txt ~{bam}
    >>>

    output {
        String model = read_string("model.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.11"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.11"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.11"
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

task Filter {
    input {
        File bam

        Int num_cpus = 2
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow filter -o ~{prefix}.filtered.bam ~{bam}
    >>>

    output {
        File filtered_bam = "~{prefix}.filtered.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.11"
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

        longbow segment ~{bam} | longbow extract -o ~{prefix}.extracted.bam
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.11"
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

task Correct {
    input {
        File bam
        String model

        Int num_cpus = 2
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    Map[String, String] allowlists = {
        "mas10v2":          "/longbow/resources/barcodes/cellranger/737K-august-2016.txt",
        "mas15v2":          "/longbow/resources/barcodes/cellranger/737K-august-2016.txt",
        "mas10threeP":      "/longbow/resources/barcodes/cellranger/3M-february-2018.txt.gz",
        "mas15threeP":      "/longbow/resources/barcodes/cellranger/3M-february-2018.txt.gz",
        "mas15teloprimev2": "/longbow/resources/barcodes/indexes/tpv2.indices.txt"
    }

    command <<<
        set -euxo pipefail

        longbow correct-tag -b CR -c CB -a ~{allowlists[model]} -o ~{prefix}.extracted.bam ~{bam}
    >>>

    output {
        File corrected_bam = "~{prefix}.corrected.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.17"
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
        File summary = "~{prefix}_summary_stats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.12"
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
