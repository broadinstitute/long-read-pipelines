version 1.0

# TODO: describe purpose
task Minimap2 {
    input {
        File shard
        File ref_fasta
        String SM
        String ID
        String PL
        Boolean? reads_are_corrected
        String docker
    }

    Boolean correct = select_first([reads_are_corrected, false])
    String map_arg = if (PL == "ONT") then "map-ont" else "map-pb"
    String correction_arg = if (correct) then "asm20" else map_arg

    Int cpus = 4
    Int disk_size = ceil(size(ref_fasta, "GB")) + 4*ceil(size(shard, "GB"))

    String aligned_shard_name = basename(shard, ".bam") + ".aligned.bam"

    command <<<
        set -euxo pipefail

        RG=`python /usr/local/bin/merge_read_group_tags.py --ID ~{ID} --SM ~{SM} --PL ~{PL} ~{shard}`
        samtools fastq ~{shard} | minimap2 -ayY --MD --eqx -x ~{correction_arg} -R ${RG} -t ~{cpus} ~{ref_fasta} - | samtools view -b - > temp.aligned.unsorted.bam
        java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RepairLongReadBam -I ~{shard} -A temp.aligned.unsorted.bam -O temp.aligned.unsorted.repaired.bam -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
        samtools sort -@~{cpus} -m4G -o ~{aligned_shard_name} temp.aligned.unsorted.repaired.bam
    >>>

    output {
        File aligned_shard = "~{aligned_shard_name}"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}