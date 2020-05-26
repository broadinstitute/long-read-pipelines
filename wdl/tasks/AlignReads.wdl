version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.7/wdl/tasks/Structs.wdl"

task Minimap2 {
    input {
        Array[File] reads
        File ref_fasta

        String RG
        String map_preset

        String prefix = "out"
        Int cpus = 4
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        MAP_PARAMS="-ayYL --MD --eqx -x ~{map_preset} -R ~{RG} -t ~{cpus} ~{ref_fasta}"
        FILE="~{reads[0]}"
        FILES="~{sep=' ' reads}"

        # We write to a SAM file before sorting and indexing because rarely, doing everything
        # in a single one-liner leads to a truncated file error and segmentation fault of unknown
        # origin.  Separating these commands requires more resources, but is more reliable overall.

        if [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            cat $FILES | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            zcat $FILES | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta$ ]] || [[ "$FILE" =~ \.fa$ ]]; then
            cat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta.gz$ ]] || [[ "$FILE" =~ \.fa.gz$ ]]; then
            zcat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.bam$ ]]; then
            samtools fastq $FILES | minimap2 $MAP_PARAMS - > tmp.sam
        else
            echo "Did not understand file format for '$FILE'"
            exit 1
        fi

        samtools sort -@~{cpus} -m4G --no-PG -o ~{prefix}.bam tmp.sam
        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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

task MergeBams {
    input {
        Array[File] bams

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bams, "GB"))

    command <<<
        set -euxo pipefail

        samtools merge -p -c -@2 --no-PG ~{prefix}.bam ~{sep=" " bams}
        samtools index ~{prefix}.bam
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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
