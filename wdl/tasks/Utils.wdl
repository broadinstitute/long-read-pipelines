version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test-2/wdl/tasks/Structs.wdl"

task ShardLongReads {
    input {
        Array[String] unmapped_files
        Int? num_reads_per_split

        RuntimeAttr? runtime_attr_override
    }

    Int nr = select_first([num_reads_per_split, 200000])
    Int disk_size = 4*ceil(size(unmapped_files, "GB"))
    Int num_files = length(unmapped_files)

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -jar /usr/local/bin/gatk.jar ShardLongReads -I ~{sep=' -I ' unmapped_files} -nr ~{nr} -O ./ -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        Array[File] unmapped_shards = glob("*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-utils:0.01.05"
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

task PrepareManifest {
    input {
        Array[String] files

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        echo ~{sep=' ' files} | sed 's/ /\n/g' > manifest.txt
    >>>

    output {
        File manifest = "manifest.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-utils:0.01.05"
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

task EchoManifest {
    input {
        File manifest

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        cat ~{manifest}
    >>>

    output {
        File out = stdout()
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-utils:0.01.05"
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

task ChunkManifest {
    input {
        File manifest
        Int manifest_lines_per_chunk

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        split -a 5 -d --additional-suffix=".txt" -l ~{manifest_lines_per_chunk} ~{manifest} chunk_
    >>>

    output {
        Array[File] manifest_chunks = glob("chunk_*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-utils:0.01.05"
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

# TODO: investigate if samtools is better for this task
# Sort BAM file by coordinate order
# copied from "dsde_pipelines_tasks/BamProcessing.wdl", with
# customization on the runtime block, and "preemptible_tries" taken out
task SortSam {
    input {
        File input_bam
        String output_bam_basename
        Int compression_level

        RuntimeAttr? runtime_attr_override
    }

    # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
    # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
    Float sort_sam_disk_multiplier = 3.25
    Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

    command {
        java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
            SortSam \
            INPUT=~{input_bam} \
            OUTPUT=~{output_bam_basename}.bam \
            SORT_ORDER="coordinate" \
            CREATE_INDEX=true \
            CREATE_MD5_FILE=true \
            MAX_RECORDS_IN_RAM=300000 \
            VALIDATION_STRINGENCY=SILENT
    }

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bam_index = "~{output_bam_basename}.bai"
        File output_bam_md5 = "~{output_bam_basename}.bam.md5"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             5,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
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

task MakeChrIntervalList {
    input {
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | awk '{ print $2 "\t" 1 "\t" $3 }' | sed 's/[SL]N://g' | grep -v -e random -e chrUn -e decoy -e alt -e HLA -e EBV > chrs.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("chrs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-metrics:0.01.07"
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

task FastaToSam {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Float fasta_sam_disk_multiplier = 3.25
    Int disk_size = ceil(fasta_sam_disk_multiplier * size(fasta, "GiB")) + 20

    command <<<
        python /usr/local/bin/prepare_run.py ~{fasta}
    >>>

    output {
        File output_bam = "unmapped.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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

task CountFastqRecords {
    input {
        File fastq

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(2 * size(fastq, "GiB"))

    command <<<
        set -euxo pipefail

        FILE="~{fastq}"
        if [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            cat ~{fastq} | awk '{s++}END{print s/4}'
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            zcat ~{fastq} | awk '{s++}END{print s/4}'
        fi
    >>>

    output {
        Int num_records = read_int(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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

task CountFastaRecords {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(2 * size(fasta, "GiB"))

    command <<<
        grep -c '>' ~{fasta}
    >>>

    output {
        Int num_records = read_int(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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

task CountBamRecords {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(2 * size(bam, "GiB"))

    command <<<
        samtools view ~{bam} | wc -l
    >>>

    output {
        Int num_records = read_int(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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

task GrepCountBamRecords {
    input {
        String bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-vc" else "-c"

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ~{samfilter} ~{bam} | grep ~{arg} ~{regex} > ~{prefix}.txt
    >>>

    output {
        Int num_records = read_int("~{prefix}.txt")
        File num_records_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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

task GrepCountUniqueBamRecords {
    input {
        String bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-v" else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ~{samfilter} ~{bam} | grep ~{arg} ~{regex} | > ~{prefix}.txt
    >>>

    output {
        Int num_records = read_int("~{prefix}.txt")
        File num_records_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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

task Sum {
    input {
        Array[Int] ints
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python -c "print(~{sep="+" ints})" > ~{prefix}.txt
    >>>

    output {
        Int sum = read_int("~{prefix}.txt")
        File sum_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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

task BamToTable {
    input {
        String bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ~{bam} | perl -n -e '($nm) = $_ =~ /NM:i:(\d+)/; ($as) = $_ =~ /AS:i:(\d+)/; ($za) = $_ =~ /ZA:Z:(\w+|\.)/; ($zu) = $_ =~ /ZU:Z:(\w+|\.)/; ($cr) = $_ =~ /CR:Z:(\w+|\.)/; ($cb) = $_ =~ /CB:Z:(\w+|\.)/; @a = split(/\s+/); print join("\t", $a[0], $a[1], $a[2], $a[3], $a[4], length($a[9]), $nm, $as, $za, $zu, $cr, $cb, $a[1], ($a[1] & 0x1 ? "paired" : "unpaired"), ($a[1] & 0x4 ? "unmapped" : "mapped"), ($a[1] & 0x10 ? "rev" : "fwd"), ($a[1] & 0x100 ? "secondary" : "primary"), ($a[1] & 0x800 ? "supplementary" : "non_supplementary")) . "\n"' | gzip > ~{prefix}.txt.gz
    >>>

    output {
        File table = "~{prefix}.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/broad-long-read-pipelines/lr-align:0.01.25"
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
