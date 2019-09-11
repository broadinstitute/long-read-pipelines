version 1.0

import "Structs.wdl"

# TODO: describe purpose
task Minimap2 {
    input {
        File shard
        File ref_fasta
        String SM
        String ID
        String PL
        Boolean? reads_are_corrected
        
        RuntimeAttr? runtime_attr_override
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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          "~{cpus}", 
        mem_gb:             20, 
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "kgarimella/lr-align:0.01.17"
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