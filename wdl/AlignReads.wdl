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
        Boolean? reads_are_rna

        RuntimeAttr? runtime_attr_override
    }

    Boolean correct = select_first([reads_are_corrected, false])
    Boolean rna_reads = select_first([reads_are_rna, false])
    String map_arg = if (PL == "ONT") then "map-ont" else "map-pb"
    String correction_arg = if (correct) then "asm20" else map_arg
    String map_preset = if (rna_reads) then "splice" else correction_arg

    Int cpus = 4
    Int disk_size = ceil(size(ref_fasta, "GB")) + 10*ceil(size(shard, "GB"))

    String aligned_shard_name = basename(shard, ".bam") + ".aligned.bam"

    command <<<
        set -euxo pipefail

        RG=`python /usr/local/bin/merge_read_group_tags.py --ID ~{ID} --SM ~{SM} --PL ~{PL} ~{shard}`
        samtools fastq ~{shard} | minimap2 -ayYL --cs --MD --eqx -x ~{map_preset} -R ${RG} -t ~{cpus} ~{ref_fasta} - | samtools view -b - > temp.aligned.unsorted.bam
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
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.18"
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

# A simple task to covert SAM-format alignment to PAF format
task SAMtoPAF {
    input {
        File sam_formatted_file
        File? index

        RuntimeAttr? runtime_attr_override
    }

    # currently we only support bam (not sam or cram)
    String prefix = basename(sam_formatted_file, ".bam")

    Int disk_size = 2*ceil(size(sam_formatted_file, "GB"))

    command <<<
        set -eu

        filename=$(basename -- ~{sam_formatted_file})
        extension="${filename##*.}"
        if [[ "$extension" == "sam" ]]; then
            /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                ~{sam_formatted_file} \
                > ~{prefix}".paf"
        elif [[ "$extension" == "bam" ]]; then
            samtools view -h ~{sam_formatted_file} | \
                /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                - \
                > ~{prefix}".paf"
        else
            echo "Currently we only support SAM or BAM (not CRAM)." && exit 1;
        fi
    >>>

    output {
        File pat_formatted_file = "~{prefix}.paf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.18"
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