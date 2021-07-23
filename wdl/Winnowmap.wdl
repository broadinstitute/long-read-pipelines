version 1.0

######################################################################################
## A pipeline for running the Winnowmap aligner
######################################################################################

import "tasks/Structs.wdl"

workflow Winnowmap {
    input {
        File reads
        String prefix
        String experiment_type
        File ref_fasta
        File ref_meryl_kmers
    }

    parameter_meta {
        reads:           "Raw reads in either fa or fq format"
        prefix:          "Output file prefix"
        experiment_type: "CLR, CCS, ONT, asm5, asm10, or asm20"
    }

    Map[String, String] map_presets = {
        'CLR':    'map-pb-clr',
        'CCS':    'map-pb',
        'ONT':    'map-ont',
        'asm5':   'asm5',
        'asm10':   'asm10',
        'asm20':   'asm20'
    }

    String preset = map_presets[experiment_type]

    call Align { input: reads = reads, prefix = prefix, preset = preset, ref_fasta = ref_fasta, ref_meryl_kmers = ref_meryl_kmers }

    output {
        File aligned_bam = Align.aligned_bam
        File aligned_bai = Align.aligned_bai
    }
}

task Align {
    input {
        File reads
        String prefix
        String preset
        File ref_fasta
        File ref_meryl_kmers

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(reads, "GB")) + ceil(size(ref_fasta, "GB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        winnowmap -W ~{ref_meryl_kmers} -t $num_core -ax ~{preset} ~{ref_fasta} ~{reads} | samtools sort -@ $num_core -O BAM -o ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-winnowmap:2.03"
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