version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.28/wdl/tasks/Structs.wdl"

task IndexUnalignedBam {
    input {
        String input_bam

        RuntimeAttr? runtime_attr_override
    }

    String bri_path = basename(input_bam) + ".bri"

    command <<<
        set -euxo pipefail

        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        HTS_AUTH_LOCATION=/tmp/token_fifo bri index -v -i ~{bri_path} ~{input_bam}
    >>>

    output {
        File bri = bri_path
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.22"
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

task MakeReadNameManifests {
    input {
        File input_bri
        Int N = 300

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        bri show ~{input_bri} | awk -F'/' 'BEGIN { getline; zmw=$2; line=$0 } { if ($2 != zmw) { print line; line = $0; } else { line = line "," $0; } zmw=$2; } END { print line; }' | sort -t'/' -n -k2 > read_list.txt
        split -a 5 -d --additional-suffix=".txt" -n l/~{N} read_list.txt chunk_
    >>>

    output {
        File manifest_full = "read_list.txt"
        Array[File] manifest_chunks = glob("chunk_*.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.22"
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


task ExtractReadsInManifest {
    input {
        String input_bam
        File input_bri
        File read_name_manifest

        Int batch_size = 10000

        RuntimeAttr? runtime_attr_override
    }

    # Estimate disk size needed assuming each read is ~100 KB (a deliberate overestimate)
    Int typical_read_disk_size = 100
    Int num_reads = length(read_lines(read_name_manifest))
    Int disk_size = ceil(2*typical_read_disk_size*num_reads/1000000)

    command <<<
        set -ux

        # We'll batch our fetches into num_reads/batch_size requests.
        sed 's/,/\n/g' ~{read_name_manifest} | sort -t'/' -n -k2 > readnames.txt
        split -a 5 -d --additional-suffix=".txt" -l ~{batch_size} readnames.txt subchunk_

        # Get renewable auth token; see comment from @daviesrob at https://github.com/samtools/htslib/issues/803 .
        mkfifo /tmp/token_fifo
        ( while true ; do curl --retry 3 -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &

        # Fetch read batches in parallel, staggering requests a bit so that the Google access token website doesn't get
        # hammered with requests.  We'll also abort and restart if fetching a small batch of reads takes too long.
        HTS_AUTH_LOCATION=/tmp/token_fifo samtools view --no-PG -H ~{input_bam} > header.sam
        parallel --delay 30 --retries 3 --timeout 1800 'test $(cat {} | wc -l) -le $(cat {} | HTS_AUTH_LOCATION=/tmp/token_fifo bri get -i ~{input_bri} ~{input_bam} - | samtools view --no-PG -b > {.}.bam && samtools view {.}.bam | wc -l)' ::: subchunk_*

        # Merge all the pieces together.
        samtools cat --no-PG -h header.sam -o reads.bam subchunk_*.bam

        # Because we may process aligned data as well (where a read can appear more than once), we determine success
        # to be when the actual number of reads we extracted is equal to or greater than what we expected, rather
        # than strictly requiring each to be equal to one another.
        EXP_NUM_READS=$(cat readnames.txt | wc -l)
        ACT_NUM_READS=$(samtools view reads.bam | wc -l)
        echo $EXP_NUM_READS $ACT_NUM_READS

        if [ "$EXP_NUM_READS" -gt "$ACT_NUM_READS" ]
        then
            rm reads.bam
            exit 1
        fi

        exit 0
    >>>

    output {
        File reads = "reads.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  5,
        max_retries:        3,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.22"
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
