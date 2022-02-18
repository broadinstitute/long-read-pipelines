version 1.0

import "Structs.wdl"

workflow Hifiasm {
    input {
        File reads
        String prefix
        Boolean trio = false
        File? mat_yak
        File? pat_yak
    }

    if (trio) {
        call Assemble_trio as Asm {
            input:
                reads  = reads,
                prefix = prefix,
                mat_yak = mat_yak,
                pat_yak = pat_yak
        }
    }

    if (!trio) {
        call Assemble as Asm {
            input:
                reads  = reads,
                prefix = prefix,
        }
    }
}

task Assemble_trio {
    input {
        File reads
        String prefix = "out"
        File mat_yak
        File pat_yak
        Int num_cpus = 32

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
        num_cpus: "number of CPUs to parallelize over"
        mat_yak:  "maternal yak file"
        pat_yak:  "paternal yak file"
    }

    Int disk_size = 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        hifiasm -o ~{prefix} -t~{num_cpus} -1 ~{pat_yak} -2 ~{mat_yak} ~{reads}
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.dip.hap1.p_ctg.gfa > ~{prefix}.dip.hap1.p_ctg.fa
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.dip.hap2.p_ctg.gfa > ~{prefix}.dip.hap2.p_ctg.fa
    >>>

    output {
        File h1_fa = "~{prefix}.dip.hap1.p_ctg.fa"
        File h2_fa = "~{prefix}.dip.hap2.p_ctg.fa"
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

        Int num_cpus = 32

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
        num_cpus: "number of CPUs to parallelize over"
    }

    Int disk_size = 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        hifiasm -o ~{prefix} -t~{num_cpus} ~{reads}
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.p_ctg.gfa > ~{prefix}.p_ctg.fa
    >>>

    output {
        File gfa = "~{prefix}.p_ctg.gfa"
        File fa = "~{prefix}.p_ctg.fa"
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
