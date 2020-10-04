version 1.0

task RaconPolish {
    input {
        File reads
        File draft_assembly

        Int n_rounds
    }

    Int mem_size = 4 * ceil(size(reads, "GB") + size(draft_assembly, "GB"))
    Int disk_size = mem_size

    command <<<
        set -euxo pipefail

        cp ~{draft_assembly} input_draft.fasta
        for i in {1..~{n_rounds}}
        do
          minimap2 -ax map-ont input_draft.fasta ~{reads} > aln.sam
          racon -c 1 -m 8 -x -6 -g -8 -w 500 -t 8 ~{reads} aln.sam input_draft.fasta > polished_${i}_draft.fasta
          cp polished_${i}_draft.fasta input_draft.fasta
        done
    >>>

    output {
        File final_polished_assembly = "input_draft.fasta"
        Array[File] incremental_polished_assemblies = glob("polished_*_draft.fasta")
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
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-east1-c"]
        cpuPlatform:            "Intel Haswell"
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-racon:0.1.0"
    }

}
