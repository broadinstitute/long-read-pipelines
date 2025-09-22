version 1.0

import "../../structs/Structs.wdl"

workflow Hifiasm {

    meta {
        description: "We run two HiFiasm jobs, one for getting alternative contigs and one for getting the haplotigs. And we take the primary assembly from the first job."
    }
    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
    }

    input {
        File reads
        String prefix

        Array[String] zones = ["us-central1-c", "us-central1-f", "us-central1-a", "us-central1-b"]
    }

    call AssembleForAltContigs {
        input:
            reads  = reads,
            prefix = prefix,
            zones = zones
    }

    call AssembleForHaplotigs {
        input:
            reads  = reads,
            prefix = prefix,
            zones = zones
    }

    output {
        File primary_gfa  = AssembleForAltContigs.primary_gfa
        File primary_tigs = AssembleForAltContigs.primary_tigs

        File alternate_gfa  = AssembleForAltContigs.alternate_gfa
        File alternate_tigs = AssembleForAltContigs.alternate_tigs

        File log_in_pVSa_mode = AssembleForAltContigs.log

        ###########
        Array[File] phased_gfas = AssembleForHaplotigs.phased_gfas
        Array[File] phased_tigs = AssembleForHaplotigs.phased_tigs

        File log_in_hap_mode = AssembleForHaplotigs.log

        # these two are saved, but the one generated in the primary VS alternate mode are preferred
        File primary_gfa_in_hap_mode  = AssembleForHaplotigs.primary_gfa
        File primary_tigs_in_hap_mode = AssembleForHaplotigs.primary_fa
    }
}

task AssembleForHaplotigs {
    input {
        File reads
        String prefix = "out"
        Array[String] zones

        RuntimeAttr? runtime_attr_override
    }

    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else proposed_memory  # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int disk_size = 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        time hifiasm \
            -o ~{prefix} \
            -t~{num_cpus} \
            ~{reads} \
            2>&1 | tee hifiasm.log

        tree -h .

        # GFA graph to contigs, primary
        # outputs generated this way has "bp" in their names
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.bp.p_ctg.gfa \
            > ~{prefix}.bp.p_ctg.fa

        ls "~{prefix}.bp.hap"*".p_ctg.gfa"

        # GFA graph to contigs, for each haplotig set
        for haplotype_gfa in ~{prefix}.bp.hap*.p_ctg.gfa; do
            filename=$(basename -- "${haplotype_gfa}")
            haplotype="${filename%.*}"
            awk '/^S/{print ">"$2; print $3}' \
                "${haplotype_gfa}" \
                > "${haplotype}".fa
        done
    >>>

    output {
        # these are saved, but the one with alt contigs genearted will be preferred for now
        File primary_gfa = "~{prefix}.bp.p_ctg.gfa"
        File primary_fa = "~{prefix}.bp.p_ctg.fa"

        Array[File] phased_gfas = glob("~{prefix}.bp.hap*.p_ctg.gfa")
        Array[File] phased_tigs = glob("~{prefix}.bp.hap*.p_ctg.fa")

        File log = "hifiasm.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.16.1"
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
        zones:                  "${sep=' ' zones}"
    }
}

task AssembleForAltContigs {
    input {
        File reads
        String prefix = "out"
        Array[String] zones

        RuntimeAttr? runtime_attr_override
    }
    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else if proposed_memory > 512 then 512 else proposed_memory # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int disk_size = 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        time hifiasm \
            -o ~{prefix} \
            -t~{num_cpus} \
            --primary \
            ~{reads} \
            2>&1 | tee hifiasm.log

        tree -h .

        # tricky, outputs generated this way has no "bp" in their file names
        # GFA graph to contigs, primary
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.p_ctg.gfa \
            > ~{prefix}.p_ctg.fa

        # GFA graph to contigs, alternate
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.a_ctg.gfa \
            > ~{prefix}.a_ctg.fa
    >>>

    output {
        File primary_gfa  = "~{prefix}.p_ctg.gfa"
        File primary_tigs = "~{prefix}.p_ctg.fa"

        File alternate_gfa  = "~{prefix}.a_ctg.gfa"
        File alternate_tigs = "~{prefix}.a_ctg.fa"

        File log = "hifiasm.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.16.1"
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
        zones:                  "${sep=' ' zones}"
    }
}
