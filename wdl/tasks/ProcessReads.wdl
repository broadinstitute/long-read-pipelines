version 1.0

import "Structs.wdl"

workflow ProcessReads {
    input {
        File manifest_chunk
        File ref_fasta
        Map[String, String] run_info
        Boolean? reads_are_rna
    }

    call PreprocessReads {
        input:
            manifest = manifest_chunk,
            platform = run_info['PL']
    }

    call Minimap2 as AlignReads {
        input:
            shard = PreprocessReads.unmapped_shard,
            ref_fasta = ref_fasta,
            SM = run_info['SM'],
            ID = run_info['ID'],
            PL = run_info['PL'],
            reads_are_corrected = true,
            reads_are_rna = reads_are_rna
    }

    output {
        File remaining_shard = PreprocessReads.remaining_shard
        File aligned_shard = AlignReads.aligned_shard
    }
}

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

    Int cpus = 8
    Int disk_size = 100

    String aligned_shard_name = basename(shard, ".bam") + ".aligned.bam"

    command <<<
        set -euxo pipefail

        df -h .

        RG=`python /usr/local/bin/merge_read_group_tags.py --ID ~{ID} --SM ~{SM} --PL ~{PL} ~{shard}`
        samtools fastq ~{shard} | minimap2 -ayY --MD --eqx -x ~{map_preset} -R ${RG} -t ~{cpus} ~{ref_fasta} - | samtools view -b - > temp.aligned.unsorted.bam

        df -h .

        java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RepairLongReadBam -I ~{shard} -A temp.aligned.unsorted.bam -O temp.aligned.unsorted.repaired.bam -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
        samtools sort -@~{cpus} -m4G -o ~{aligned_shard_name} temp.aligned.unsorted.repaired.bam

        df -h .
    >>>

    output {
        File aligned_shard = "~{aligned_shard_name}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.21"
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

task Minimap2Fastq {
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

    Int cpus = 8
    Int disk_size = 100

    String aligned_shard_name = basename(shard, ".bam") + ".aligned.bam"

    command <<<
        set -euxo pipefail

        cat ~{shard} | python3 /usr/local/bin/cat_as_fastq.py | minimap2 -ayY --MD --eqx -x ~{map_preset} -R '@RG\tID:~{ID}\tSM:~{SM}\tPL:~{PL}' -t ~{cpus} ~{ref_fasta} - | samtools sort -@~{cpus} -m4G -o ~{aligned_shard_name} -
    >>>

    output {
        File aligned_shard = "~{aligned_shard_name}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.22"
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

task PreprocessReads {
    input {
        File manifest
        String platform

        RuntimeAttr? runtime_attr_override
    }

    Array[File] files = read_lines(manifest)

    Int max_length = 21000
    Int min_passes = 2
    Int cpus = 4
    Int disk_size = 5*ceil(size(files[0], "GB")*length(files))

    command <<<
        set -euxo pipefail

        PLATFORM="~{platform}"

        if [ $PLATFORM == "ONT" ]
        then
            python /usr/local/bin/prepare_run.py ~{sep=' ' files}
            samtools view -H unmapped.bam | samtools view -b > remaining.bam

            echo "ZMWs input          (A)  : 0
ZMWs generating CCS (B)  : 0 (100.00%)
ZMWs filtered       (C)  : 0 (0.00%)

Exclusive ZMW counts for (C):
No usable subreads       : 0 (0.00%)
Below SNR threshold      : 0 (0.00%)
Lacking full passes      : 0 (0.00%)
Heteroduplexes           : 0 (0.00%)
Min coverage violation   : 0 (0.00%)
Draft generation error   : 0 (0.00%)
Draft above --max-length : 0 (0.00%)
Draft below --min-length : 0 (0.00%)
Lacking usable subreads  : 0 (0.00%)
CCS did not converge     : 0 (0.00%)
CCS below minimum RQ     : 0 (0.00%)
Unknown error            : 0 (0.00%)" > ccs_report.txt
        else
            ccs --max-length ~{max_length} --min-passes ~{min_passes} -j ~{cpus} ~{files[0]} unmapped.bam

            java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RecoverUncorrectedReads -I ~{files[0]} -C unmapped.bam -O remaining.bam -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
        fi
    >>>

    output {
        File unmapped_shard = "unmapped.bam"
        File remaining_shard = "remaining.bam"
        File report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             40,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.21"
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
        Array[File] shards
        String merged_name

        RuntimeAttr? runtime_attr_override
    }

    String merged_bai = basename(merged_name, ".bam") + ".bai"
    Int disk_size = 3*ceil(size(shards, "GB"))

    command <<<
        set -euxo pipefail

        samtools merge -@2 ~{merged_name} ~{sep=" " shards}
        samtools index -@2 ~{merged_name} ~{merged_bai}
    >>>

    output {
        File merged = "~{merged_name}"
        File merged_bai = "~{merged_bai}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.21"
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

task FastqsToUnmappedBam {
    input {
        Array[File] manifest_chunks
        String ID
        String SM
        String PL

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1000

    command <<<
        set -euxo pipefail

        mkdir tmp
        cat ~{sep=' ' manifest_chunks} | gsutil -m cp -I tmp/
        ls -lah tmp/

        echo '@HD VN:1.5 SO:unknown' | sed 's/ /\t/g' > header.sam
        echo '@RG ID:~{ID} SM:~{SM} PL:~{PL}' | sed 's/ /\t/g' >> header.sam

        ((cat header.sam) && (cat tmp/*.fastq | paste - - - - | sed "s/^@//" | awk -F "\\t" '{ sub(/ .+/, "", $1); print $1, "4", "*", "0", "0", "*", "*", "0", "0", $2, $4, "RG:Z:~{ID}" }' | sed "s/ /\\t/g")) | samtools view -b > unmapped.bam
    >>>

    output {
        File unmapped = "unmapped.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-utils:0.01.02"
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

task ValidateBam {
    input {
        File input_bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(1.2*size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx4g -jar /usr/local/bin/gatk.jar ValidateSamFile -I ~{input_bam} -O bam_validation_report.txt --IGNORE_WARNINGS
    >>>

    output {
        File report = "bam_validation_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.21"
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

