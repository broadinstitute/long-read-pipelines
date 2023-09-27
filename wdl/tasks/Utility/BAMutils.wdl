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

task FilterBamByLen {
    meta {
        desciption: "Filter a BAM by sequence length"
        note: "It's assumed that for aligned BAM, aligned was done without hard clipping turned on. If this assumption isn't met, the resulting BAM may be corrupt."
    }
    parameter_meta {
        len_threshold_inclusive: "Reads longer than or equal to this length will be included."
        bam : {
            localization_optional: true
        }
    }
    input {
        File bam
        File? bai
        Int len_threshold_inclusive

        Boolean compute_yield = false
        String disk_type = "HDD"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + 2 * ceil(size(bam, "GiB"))

    String base = basename(bam, ".bam")
    String out_prefx = "~{base}.RL_ge_~{len_threshold_inclusive}"

    Boolean is_aligned = defined(bai)

    String local_bam = "/cromwell_root/~{base}.bam"

    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}
        if ~{defined(bai)}; then mv ~{bai} "~{local_bam}.bai"; touch "~{local_bam}.bai"; fi

        # get total yield in the background
        if ~{compute_yield}; then
            samtools view -@1 \
                ~{true='-F2304' false=' ' is_aligned} \
                ~{local_bam} \
            | awk -F '\t' '{print length($10)}' \
            > all.read.lengths.txt &
        fi

        if ~{is_aligned} ; then
            # note here that 2048 and 256 reads are no longer than the primary record,
            # so if the primary is already shorter than the threshold, they should be excluded too
            # on the other hand, it is possible the 2048 records are shorter that the actual query when hardclipping is turned on
            # hence we requre the input doesn't have hardclipping turned on
            samtools view -@1 -h \
                --write-index \
                -e "length(seq)>=~{len_threshold_inclusive}" \
                -o "~{out_prefx}.bam##idx##~{out_prefx}.bam.bai" \
                ~{local_bam}
        else
            samtools view -@1 -h \
                -e "length(seq)>=~{len_threshold_inclusive}" \
                -o "~{out_prefx}.bam" \
                ~{local_bam}
        fi

        if ~{compute_yield}; then
            samtools view -@1 \
                ~{true='-F2304' false=' ' is_aligned} \
                "~{out_prefx}.bam" \
            | awk -F '\t' '{print length($10)}' \
            > filtered.read.lengths.txt

            # see https://stackoverflow.com/questions/450799/shell-command-to-sum-integers-one-per-line
            awk '{s+=$1} END {printf "%.0f", s}' filtered.read.lengths.txt \
            > filtered.yield.txt

            wait # make sure total yield gather in the background is done
            awk '{s+=$1} END {printf "%.0f", s}' all.read.lengths.txt \
            > all.yield.txt
        fi
    >>>
    output {
        File  fBAM = "~{out_prefx}.bam"
        File? fBAI = "~{out_prefx}.bam.bai"

        Float?    total_yield = read_float("all.yield.txt")
        Float? filtered_yield = read_float("filtered.yield.txt")
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             4,
        mem_gb:                16,
        disk_gb:               disk_size,
        preemptible_tries:     1,
        max_retries:           1,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task InferSampleName {
    meta {
        description: "Infer sample name encoded on the @RG line of the header section. Fails if multiple values found, or if SM ~= unnamedsample."
    }

    input {
        File bam
        File? bai
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

task BamToFastq {

    meta {
        description : "Convert a BAM file to a fastq file."
    }

    parameter_meta {
        bam: {localization_optional: true}
        prefix: "Prefix for the output fastq file."
        disk_type: "type of disk to use"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        String prefix

        String disk_type = "HDD"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 3 * ceil(size(bam, "GiB"))

    String base = basename(bam)
    String local_bam = "/cromwell_root/~{base}"
    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}
        time samtools fastq -@1 -t -0 ~{prefix}.fq.gz ~{local_bam}
    >>>

    output {
        File reads_fq = "~{prefix}.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.1"
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
