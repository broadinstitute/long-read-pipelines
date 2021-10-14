version 1.0

import "Structs.wdl"

workflow Hifiasm {
    input {
        File reads
        File? mother_bam
        File? father_bam
        String prefix
    }

    if (defined(mother_bam) && defined(father_bam)) {
        call Yak as YakMother { input: bam = select_first([mother_bam]) }
        call Yak as YakFather { input: bam = select_first([father_bam]) }
    }

    call Assemble {
        input:
            reads  = reads,
            prefix = prefix,
            yak_mother = YakMother.yak,
            yak_father = YakFather.yak,
    }

    scatter (gfa in Assemble.gfa) {
        call GfaToFa { input: gfa = gfa }
    }

    output {
        File gfa = Assemble.gfa[0]
        File fa = GfaToFa.fa[0]
        Array[File] haps = GfaToFa.fa
    }
}

task Yak {
    input {
        File bam
        String prefix = "out"

        Int num_cpus = 16

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:      "reads in BAM format"
        prefix:   "prefix to apply to assembly output filenames"
        num_cpus: "number of CPUs to parallelize over"
    }

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        yak count -b37 -t16 -o ~{prefix}.yak <(samtools fastq ~{bam})
    >>>

    output {
        File yak = "~{prefix}.yak"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             150,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.13"
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

task Assemble {
    input {
        File reads
        String prefix = "out"
        File? yak_mother
        File? yak_father

        Int num_cpus = 64

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:        "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:       "prefix to apply to assembly output filenames"
        yak_mother:   "Yak file for mother's reads"
        yak_father:   "Yak file for father's reads"
        num_cpus:     "number of CPUs to parallelize over"
    }

    Int disk_size = 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        hifiasm -o ~{prefix} \
                -t~{num_cpus} \
                ~{true='-1' false='' defined(yak_mother)} ~{select_first([yak_mother, ''])} \
                ~{true='-2' false='' defined(yak_father)} ~{select_first([yak_father, ''])} \
                ~{reads}
    >>>

    output {
        Array[File] gfa = glob("*.p_ctg.gfa")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             150,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.13"
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

task GfaToFa {
    input {
        File gfa

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2 * ceil(size(gfa, "GB"))
    String prefix = basename(gfa, ".gfa")

    command <<<
        set -euxo pipefail

        awk '/^S/{print ">"$2; print $3}' ~{gfa} > ~{prefix}.fa
    >>>

    output {
        File fa = "~{prefix}.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.13"
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
