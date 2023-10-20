version 1.0

import "Structs.wdl"

task GetDefaultDir {
    input {
        String workflow_name

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        NAME=$(cat gcs_localization.sh | grep 'source bucket' | sed 's/# Localize files from source bucket //' | sed 's/ to container.*//' | sed "s/'//g")

        echo "gs://$NAME/results/~{workflow_name}"
    >>>

    output {
        String path = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task SortBam {
    input {
        File input_bam
        String prefix = "sorted"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        input_bam: "input BAM"
        prefix:    "[default-valued] prefix for output BAM"
    }

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        samtools sort -@$num_core -o ~{prefix}.bam ~{input_bam}
        samtools index ~{prefix}.bam
    >>>

    output {
        File sorted_bam = "~{prefix}.bam"
        File sorted_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task MergeFastqGzs {
    input {
        Array[File] fastq_gzs

        String prefix = "merged"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(fastq_gzs, "GB"))

    command <<<
        cat ~{sep=' ' fastq_gzs} > ~{prefix}.fastq.gz
    >>>

    output {
        File merged_fastq_gz = "~{prefix}.fastq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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


task Cat {
    input {
        Array[File] files
        Boolean has_header = false
        String out = "out.txt"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        files:      "text files to combine"
        has_header: "files have a redundant header"
        out:        "[default-valued] output filename"
    }

    Int disk_size = 4*ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        HAS_HEADER=~{true='1' false='0' has_header}

        if [ "$HAS_HEADER" = 1 ]; then
            ((head -1 ~{files[0]}) && (cat ~{sep=' ' files} | xargs -n 1 tail -n +2)) > ~{out}
        else
            cat ~{sep=' ' files} > ~{out}
        fi
    >>>

    output {
        File combined = out
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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
        Array[String] filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | \
            awk '{ print $2 "\t" 1 "\t" $3 }' | \
            sed 's/[SL]N://g' | \
            grep -v -e '^@HD' ~{true='-e' false='' length(filter) > 0} ~{sep=" -e " filter} | \
            tee chrs.txt
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

task CountFastqRecords {
    input {
        File fastq

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(fastq, "GiB"))

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

task CountFastaRecords {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(fasta, "GiB"))

    command <<<
        grep -c '>' ~{fasta}

        exit 0
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

task CountBamRecords {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta { bam: { localization_optional: true } }

    Int disk_size = 100

    command <<<
        set -eux
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -c ~{bam} > "count.txt" 2>"error.log"
        if [[ -f "error.log" ]]; then
            if [[ -s "error.log" ]]; then echo "samtools has warn/error" && cat "error.log" && exit 1; fi
        fi
    >>>

    output {
        File? samools_error = "error.log"
        Int num_records = read_int("count.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

task FilterListOfStrings {
    meta {
        description : "Filter a given list of files by a query term (essentially pipes the query into grep)."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        Array[String] list_to_filter
        String query

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        list_to_filter: "Array of strings to filter by the query."
        query: "Term to use to filter the given strings."
    }


    command <<<
        set -euxo pipefail

        \grep "~{query}" ~{write_lines(list_to_filter)} > filtered_list.txt
    >>>

    output {
        Array[String] filtered_list = read_lines("filtered_list.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "ubuntu:hirsute-20210825"
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

task FilterReadsBySamFlags {
    meta {
        description : "Filter reads based on sam flags.  Reads with ANY of the given flags will be removed from the given dataset."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam
        String sam_flags

        String extra_args = ""

        String prefix = "filtered_reads"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:   "BAM file to be filtered."
        sam_flags: "Flags for which to remove reads.  Reads with ANY of the given flags will be removed from the given dataset."
        prefix : "[Optional] Prefix string to name the output file (Default: filtered_reads)."
    }

    Int disk_size = 20 + ceil(11 * size(bam, "GiB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools view -h -b -F ~{sam_flags} -@$np ~{extra_args} ~{bam} > ~{prefix}.bam
        samtools index -@$np ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
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

task DownsampleSam {
    meta {
        description : "Downsample the given bam / sam file using Picard/GATK's DownsampleSam tool."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam

        Float probability = 0.01
        String strategy = "HighAccuracy"
        String prefix = "downsampled_reads"

        Int random_seed = 1

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:   "BAM file to be filtered."
        probability : "[Optional] Probability that a read will be emitted (Default: 0.01)."
        strategy : "[Optional] Strategy to use to downsample the given bam file (Default: HighAccuracy)."
        prefix : "[Optional] Prefix string to name the output file (Default: downsampled_reads)."
        extra_args : "[Optional] Extra arguments to pass into DownsampleSam."
    }

    Int disk_size = 20 + ceil(11 * size(bam, "GiB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        gatk DownsampleSam --VALIDATION_STRINGENCY SILENT --RANDOM_SEED ~{random_seed} -I ~{bam} -O ~{prefix}.bam -S ~{strategy} -P ~{probability} ~{extra_args}
        samtools index -@$np ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.2.0.0"
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
        File bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-vc" else "-c"

    command <<<
        set -euxo pipefail

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

task GrepCountUniqueBamRecords {
    input {
        String bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-v" else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
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

task Uniq {
    input {
        Array[String] strings

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        sort ~{write_lines(strings)} | uniq > uniq.txt
    >>>

    output {
        Array[String] unique_strings = read_lines("uniq.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task Timestamp {
    input {
        Array[String] dummy_dependencies

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        date --iso-8601=ns > timestamp.txt
    >>>

    output {
        String timestamp = read_string("timestamp.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task GetUnalignedReads {
    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        samtools view -b -f 4 ~{bam} > ~{prefix}.unaligned_reads.bam
    >>>

    output {
        File unaligned_reads_bam = "~{prefix}.unaligned_reads.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam, "GB"))

    command <<<
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

task ConvertReads {
    input {
        File reads
        String output_format
    }

    Int disk_size = 3 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        filename=~{reads}
        input_filetype=${filename##*.}
        output_filetype=~{output_format}

        if [[ ($input_filetype == "fastq" || $input_filetype == "fq") && $output_filetype == "fasta" ]]; then
            echo "Converting $input_filetype to $output_filetype"
            seqkit fq2fa $filename -o tmp.out
        elif [ $input_filetype == $output_filetype ]; then
            echo "Input filetype is the output filetype"
            mv $filename tmp.out
        else
            echo "ConvertReads does not know how to convert $input_filetype to $output_filetype"
            exit 1
        fi

        mv tmp.out converted_reads.$output_filetype
    >>>

    output {
        File converted_reads = "converted_reads.~{output_format}"
    }

    runtime {
        cpu:                    4
        memory:                 "8 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         10
        preemptible:            2
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}

task BamToBed {
    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail
        bedtools genomecov -ibam ~{bam} -bg > ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.10"
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

task BamToFastq {
    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools fastq ~{bam} | gzip > ~{prefix}.fq.gz
    >>>

    output {
        File reads_fq = "~{prefix}.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task MergeFastqs {
    input {
        Array[File] fastqs

        String prefix = "merged"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(fastqs, "GB"))

    command <<<
        mkdir tmp
        mv -t tmp ~{sep=' ' fastqs}

        cat tmp/*.fastq.gz > ~{prefix}.fastq.gz
    >>>

    output {
        File merged_fastq = "~{prefix}.fastq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

# A utility to merge several input BAMs into a single BAM.
task MergeBams {
    input {
        Array[File] bams
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bams:   "input array of BAMs to be merged"
        prefix: "[default-valued] prefix for output BAM"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

# A utility to subset a BAM to specifed loci
task SubsetBam {
    input {
        String bam
        File bai
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:    "bam to subset"
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
    }

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task SplitBam {
    input {
        File bam
        File bai
        Array[String] filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA']

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:    "bam to split"
        bai:    "index for bam file"
        filter: "contigs to ignore"
    }

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        samtools view -H ~{bam} | \
            grep '^@SQ' | \
            grep -v -e '^@HD' ~{true='-e' false='' length(filter) > 0} ~{sep=" -e " filter} | \
            awk '{ print $2 }' | \
            sed 's/SN://' |
            parallel -j+0 "samtools view -bh -o {}.bam ~{bam} {} && samtools index {}.bam"
    >>>

    output {
        Array[File] subset_bams = glob("*.bam")
        Array[File] subset_bais = glob("*.bai")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task ListBamContigs {
    input {
        String bam

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:    "input BAM from which available contigs should be listed"
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep '^@SQ' | awk '{ print $2 }' | sed 's/SN://' > chrs.txt
    >>>

    output {
        Array[String] contigs = read_lines("chrs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task ComputeGenomeLength {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        fasta:  "FASTA file"
    }

    Int disk_size = 2*ceil(size(fasta, "GB"))

    command <<<
        set -euxo pipefail

        samtools dict ~{fasta} | \
            grep '^@SQ' | \
            awk '{ print $3 }' | \
            sed 's/LN://' | \
            awk '{ sum += $1 } END { print sum }' > length.txt
    >>>

    output {
        Float length = read_float("length.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task ListFilesOfType {
    input {
        String gcs_dir
        Array[String] suffixes

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gcs_dir:  "input directory"
        suffixes: "suffix(es) for files"
    }

    Int disk_size = 1
    String in_dir = sub(gcs_dir, "/$", "")

    command <<<
        set -x

        RET=0

        while read s; do
            gsutil ls ~{in_dir}/**$s >> files.txt
        done <~{write_lines(suffixes)}

        if [[ $(wc -l files.txt) -eq 0 ]]; then
            RET=1
        fi

        exit $RET
    >>>

    output {
        Array[String] files = read_lines("files.txt")
        File manifest = "files.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task StopWorkflow {
    input {
        String reason
    }
    command <<<
        echo -e "Workflow explicitly stopped because \n  ~{reason}." && exit 1
    >>>
    runtime {
      docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task InferSampleName {
    meta {
        description: "Infer sample name encoded on the @RG line of the header section. Fails if multiple values found, or if SM ~= unnamedsample."
    }

    input {
        File bam
        File bai
    }

    parameter_meta {
        bam: {
            localization_optional: true
        }
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt
        if ! grep -q '^@RG' header.txt; then echo "No read group line found!" && exit 1; fi

        grep '^@RG' header.txt | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g' | sort | uniq > sample.names.txt
        if [[ $(wc -l sample.names.txt) -gt 1 ]]; then echo "Multiple sample names found!" && exit 1; fi
        if grep -iq "unnamedsample" sample.names.txt; then echo "Sample name found to be unnamedsample!" && exit 1; fi
    >>>

    output {
        String sample_name = read_string("sample.names.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task CheckOnSamplenames {
    meta {
        description: "Makes sure the provided sample names are all same, i.e. no mixture of sample names"
    }

    input {
        Array[String] sample_names
    }

    command <<<
        set -eux
        n_sm=$(sort ~{write_lines(sample_names)} | uniq | wc -l | awk '{print $1}')
        if [[ ${n_sm} -gt 1 ]]; then echo "Sample mixture!" && exit 1; fi
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task DeduplicateBam {
    meta {
        description: "Utility to drop (occationally happening) duplicate records in input BAM"
    }

    input {
        File aligned_bam
        File aligned_bai

        Boolean same_name_as_input = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(aligned_bam, "GB"))

    String base = basename(aligned_bam, ".bam")
    String prefix = if (same_name_as_input) then base else (base + ".dedup")

    command <<<
        echo "==========================================================="
        echo "collecting duplicate information"
        time \
            samtools view -@ 1 "~{aligned_bam}" | \
            awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' | \
            sort | uniq -d \
            > "~{aligned_bam}".duplicates.txt
        echo "==========================================================="
        echo "de-duplicating"
        time python3 /opt/remove_duplicate_ont_aln.py \
            --prefix "~{prefix}" \
            --annotations "~{aligned_bam}".duplicates.txt \
            "~{aligned_bam}"
        echo "==========================================================="
        echo "DONE"
        samtools index "~{prefix}.bam"
    >>>

    output {
        File corrected_bam = "~{prefix}.bam"
        File corrected_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
