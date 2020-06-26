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
            | python /soft_clipper.py --split-read-prefix=x${i} --clipping-threshold=~{clipping_threshold} \
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
        docker:                 "quay.io/broad-long-read-pipelines/lr-softclipper:0.5.0"
    }
}

task SplitSoftClippedReadsAssisted {
    input {
        File reads_fastq
        File reference_fasta
        Int rounds
        Int clipping_threshold
        Int ref_conflict_threshold
        File aid_reference_fasta
    }

    Int disk_size = (14 + 6 + rounds) * ceil(size(reads_fastq, "GB"))
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

          minimap2 --eqx -ax map-ont ~{aid_reference_fasta} ${input_fn} > aid_ref.sam
          minimap2 --eqx -ax map-ont ~{reference_fasta} ${input_fn} > ref.sam

          cat ref.sam | python /soft_clipper.py \
               --clipping-threshold=~{clipping_threshold} \
               --split-read-prefix=x${i} \
               --ref=aid_ref.sam \
               --ref-diff-threshold=~{ref_conflict_threshold} \
               --write-ref-conflicts-prefix=conflicts_x${i} \
            > ~{basename}_softclipped_x${i}.fastq

        done
    >>>

    output {
       Array[File] split_reads = glob("*_softclipped_*.fastq")
       Array[File] conflicting_alignments = glob("conflicts_*.bam")
       File most_split_read = "~{basename}_softclipped_x~{rounds}.fastq"
    }

    runtime {
        cpu:                    8
        memory:                 "60 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         10
        preemptible:            0
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-softclipper:0.5.0"
    }
}
