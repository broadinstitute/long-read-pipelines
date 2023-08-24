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

task MergeBamsWithSamtools {

    meta {
        description : "Merge several input BAMs into a single BAM."
    }

    parameter_meta {
        bams: {localization_optional: true}
        out_prefix: "result file will be named <out_prefix>.bam"
    }

    input {
        Array[File] bams
        String out_prefix = "out"

        String disk_type = "LOCAL"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        mkdir -p bams_dir
        time \
        gcloud storage cp ~{sep=' ' bams} /cromwell_root/bams_dir/
        ls bams_dir

        cd bams_dir && ls ./*.bam > bams.list
        time \
        samtools merge \
            -p -c --no-PG \
            -@ 3 \
            --write-index \
            -o "~{out_prefix}.bam##idx##~{out_prefix}.bam.bai" \
            -b bams.list
        mv ~{out_prefix}.bam \
           ~{out_prefix}.bam.bai \
        /cromwell_root
    >>>

    output {
        File merged_bam = "~{out_prefix}.bam"
        File merged_bai = "~{out_prefix}.bam.bai"
    }

    Int local_ssd_sz = if size(bams, "GiB") > 150 then 750 else 375
    Int pd_sz = 10 + 3*ceil(size(bams, "GiB"))
    Int disk_size = if "LOCAL" == disk_type then local_ssd_sz else pd_sz

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  if "LOCAL" == disk_type then 1 else 0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SubsetBamToLocusLocal {
    meta {
        description: "For subsetting a BAM stored in GCS, expliciting localizing the BAM"
        note: "This is intended as last resort when streaming from buckets fails"
    }

    parameter_meta {
        interval_list_file:  "a Picard-style interval list file to subset reads with"
        interval_id:         "an ID string for representing the intervals in the interval list file"
        prefix: "prefix for output bam and bai file names"
        bam: {localization_optional: true}
    }

    input {
        File bam
        File bai

        File interval_list_file
        String interval_id
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Array[String] intervals = read_lines(interval_list_file)

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    String subset_prefix = prefix + "." + interval_id

    String local_bam = "/cromwell_root/~{basename(bam)}"

    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}
        mv ~{bai} "~{local_bam}.bai" && touch "~{local_bam}.bai"

        # see man page for what '-M' means
        samtools view \
            -bhX \
            -M \
            -@ 1 \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            ~{local_bam} "~{local_bam}.bai" \
            ~{sep=" " intervals}
    >>>

    output {
        File subset_bam = "~{subset_prefix}.bam"
        File subset_bai = "~{subset_prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"  # expensive, but much faster
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
