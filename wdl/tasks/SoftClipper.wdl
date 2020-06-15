version 1.0

task SplitSoftClippedReads {
    input {
        File reads_fastq
        File reference_fasta
        Int rounds
        Int clipping_threshold
    }

    Int disk_size = (4 + rounds) * ceil(size(reads_fastq, "GB"))
    String basename = basename(reads_fastq, ".fastq")

    command <<<
        set -euxo pipefail

        for i in {1..~{rounds}}
        do
          if [[ $i -eq 1 ]];
          then
            input_fn=~{reads_fastq}
          else
            input_fn=~{basename}_softclipped_x$((i - 1)).fastq
          fi

          minimap2 --eqx -ax map-ont ~{reference_fasta} ${input_fn} \
            | samtools view -F 0x900 - \
            | python /soft_clipper.py --clipping-threshold ~{clipping_threshold} \
            > ~{basename}_softclipped_x${i}.fastq

        done
    >>>

    output {
       Array[File] split_reads = glob("*_softclipped_*.fastq")
       File most_split_read = "~{basename}_softclipped_x~{rounds}.fastq"
    }

    runtime {
        cpu:                    8
        memory:                 "32 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         10
        preemptible:            0
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-softclipper:0.1.0"
    }
}