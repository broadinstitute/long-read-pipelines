version 1.0

task RaconPolish {
    input {
        File reads
        File draft_assembly
    }

    Int mem_size = 2 * ceil(size(reads, "GB"))
    Int disk_size = 4 * ceil(size(reads, "GB") + size(draft_assembly, "GB"))

    command <<<
        set -euxo pipefail
        minimap2 -ax map-ont ~{draft_assembly} ~{reads} > aln.sam

        # -c turns on CUDA
        racon -c 1 -m 8 -x -6 -g -8 -w 500 -t 8 ~{reads} aln.sam ~{draft_assembly} > racon_polished.fasta
    >>>

    output {
        File polished_assembly = "racon_polished.fasta"
    }

    runtime {
        cpu:                    8
        # Racon has a high memory requirement. Not sure what it is exactly but you need at least
        # the size of the generated alignment file and more
        memory:                 mem_size + " GiB"
        disks:                  "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:         30
        preemptible:            0
        maxRetries:             0
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               1
        nvidiaDriverVersion:    "418.87.00"
        zones:                  ["us-east1-c"]
        cpuPlatform:            "Intel Haswell"
        docker:                 "quay.io/broad-long-read-pipelines/lr-racon:3.6.0"
    }

}