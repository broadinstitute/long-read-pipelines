version 1.0
# racon task from longread-pipeline-task Racon.wdl file
workflow Polish_Asm_by_Racon{
    meta{
        description: "a workflow that using Racon to do polish assemblies"
    }
    input{
        File inputreads
        File assembly_hap1
        File assembly_hap2
        Int n_rounds
    }
    call RaconPolish as raconhap1{input: reads=inputreads, draft_assembly=assembly_hap1, n_rounds=n_rounds}
    call RaconPolish as raconhap2{input: reads=inputreads, draft_assembly=assembly_hap2, n_rounds=n_rounds}
    output{
        File final_polished_hap1 = raconhap1.final_polished_assembly
        File final_polished_hap2 = raconhap2.final_polished_assembly
    }
}

task RaconPolish {

    meta {
        description: "Polish a draft assembly with long reads using Racon. Recommended to run a few times."
    }
    parameter_meta {
        reads:          "long reads to polish the draft assembly with"
        draft_assembly: "draft to be polished"
        n_rounds: "Number of times to run Racon"
    }

    input {
        File reads
        File draft_assembly

        Int n_rounds
    }

    Int mem_size = 4 * ceil(size(reads, "GB") + size(draft_assembly, "GB"))
    Int disk_size = mem_size

    command <<<
        set -euxo pipefail

        #cp  input_draft.fasta
        for i in {1..~{n_rounds}}
        do
          minimap2 -t 16 -ax map-hifi ~{draft_assembly} ~{reads} > aln.sam
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