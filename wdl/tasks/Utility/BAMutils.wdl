version 1.0

import "../../structs/Structs.wdl"

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information Given a single-readgroup BAM. Will fail if the information isn't present."
    }

    parameter_meta {
        uBAM: "The input BAM file."
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{uBAM} | grep "^@RG" | tr '\t' '\n' > rh_header.txt

        for attribute in ~{sep=' ' keys}; do
            value=$(grep "^${attribute}" rh_header.txt | awk -F ':' '{print $2}')
            echo -e "${attribute}\t${value}" >> "result.txt"
        done
    >>>

    output {
        Map[String, String] read_group_info = read_map("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ExtractReadsByNameRemote {
    meta {
        desciption: "Extract reads from a BAM file based on desired read names"
    }
    parameter_meta {
        bam: {
            localization_optional: true
        }
        selBam: "BAM containing selected reads"
        failed: "If this is true, the output BAM is suspected to be corrupted. Don't use and fail your overall workflow."
    }
    input {
        File  bam
        File? bai
        File qnames
        String prefix

        RuntimeAttr? runtime_attr_override
    }
    output {
        File selBam = "~{prefix}.bam"
        File? selBai = "~{prefix}.bam.bai"
        Boolean failed = read_boolean("samtools.failed.txt")
    }

    Boolean index = defined(bai)
    Int disk_size = 10 + 2*ceil(size([bam, bai], "GB"))

    command <<<
        set -eux

        # the way this works is the following:
        # 0) relying on the re-auth.sh script to export the credentials
        # 1) perform the remote sam-view subsetting in the background
        # 2) listen to the PID of the background process, while re-auth every 1200 seconds
        source /opt/re-auth.sh
        set -euxo pipefail

        echo "false" > samtools.failed.txt

        if ~{index}; then
            # see man page for what '-M' means
            samtools view \
                -bhX \
                -M \
                -@ 1 \
                --verbosity=8 \
                --write-index \
                -N ~{qnames} \
                -o "~{prefix}.bam##idx##~{prefix}.bam.bai" \
                ~{bam} ~{bai} && exit 0 || { echo "true" > samtools.failed.txt; exit 77; } &
            pid=$!
        else
            samtools view \
                -bh \
                --verbosity=8 \
                -N ~{qnames} \
                -o "~{prefix}.bam" \
                ~{bam} && exit 0 || { echo "true" > samtools.failed.txt; exit 77; } &
                pid=$!
        fi

        set +e
        count=0
        while true; do
            sleep 1200 && date && source /opt/re-auth.sh
            count=$(( count+1 ))
            if [[ ${count} -gt 6 ]]; then exit 0; fi
            if ! pgrep -x -P $pid; then exit 0; fi
        done

    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ExtractReadsByNameLocal {
    meta {
        desciption: "Extract reads from a BAM file based on desired read names"
    }
    parameter_meta {
        selBam: "BAM containing selected reads"
    }
    input {
        File  bam
        File? bai
        File qnames
        String prefix

        RuntimeAttr? runtime_attr_override
    }
    output {
        File selBam = "~{prefix}.bam"
        File? selBai = "~{prefix}.bam.bai"
        Boolean failed = read_boolean("samtools.failed.txt")
    }

    Int disk_size = 10 + 3*ceil(size([bam, bai], "GB"))

    Boolean index = defined(bai)

    command <<<
        set -eux

        if ~{index}; then
            # see man page for what '-M' means
            samtools view \
                -bhX \
                -M \
                -@ 1 \
                --verbosity=8 \
                --write-index \
                -N ~{qnames} \
                -o "~{prefix}.bam##idx##~{prefix}.bam.bai" \
                ~{bam} ~{bai}
        else
            samtools view \
                -bh \
                -@ 1 \
                --verbosity=8 \
                -N ~{qnames} \
                -o "~{prefix}.bam" \
                ~{bam}
        fi

    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
