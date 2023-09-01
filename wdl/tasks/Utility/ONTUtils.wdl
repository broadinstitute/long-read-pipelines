version 1.0

import "../../structs/Structs.wdl"

task FindSequencingSummaryFiles {

    meta {
        description: "Find sequencing summary files in an ONT basecall directory."
    }

    parameter_meta {
        gcs_input_dir: "GCS directory containing sequencing summary files."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        for summary_file in $(gsutil ls "~{indir}/**sequencing_summary*.txt*")
        do
            DIR=$(dirname $summary_file)
            echo ${DIR}

            gsutil ls "${DIR}" | grep fastq_pass && gsutil ls "${DIR}" | grep fast5_pass

            if [ $? -eq 0 ]; then
                FASTQ_COUNT=$(gsutil ls "${DIR}/fastq_pass/*.fastq*" | wc -l)
                FAST5_COUNT=$(gsutil ls "${DIR}/fast5_pass/*.fast5*" | wc -l)

                echo "${FASTQ_COUNT} ${FAST5_COUNT}"

                if [ ${FASTQ_COUNT} -eq ${FAST5_COUNT} ]; then
                    echo $summary_file >> summaries.txt
                else
                    echo "# fastq != # fast5.  Skipped."
                fi
            else
                echo "No passing fastq and fast5 files.  Skipped."
            fi

            echo ""
        done
    >>>

    output {
        Array[String] summary_files = read_lines("summaries.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
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

task GetRunInfo {

    meta{
        description: "Get ONT run info from a final summary file."
    }

    parameter_meta {
        final_summary: "Sequencing summary file."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String final_summary

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        gsutil cat "~{final_summary}" | sed 's/=[[:space:]]*$/=unknown/' | sed 's/=/\t/g' > run_info.txt
    >>>

    output {
        Map[String, String] run_info = read_map("run_info.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
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

task ListFiles {

    meta {
        description: "List files in a GCS directory."
    }

    parameter_meta {
        sequencing_summary: "Sequencing summary file."
        suffix: "Suffix of files to list."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String sequencing_summary
        String suffix

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(sub(sequencing_summary, basename(sequencing_summary), ""), "/$", "")

    command <<<
        set -euxo pipefail

        gsutil ls "~{indir}/**.~{suffix}*" | grep -v fail > files.txt
        cat files.txt | wc -l > lc.txt
    >>>

    output {
        File manifest = "files.txt"
        Int count = read_int("lc.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
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

task PartitionManifest {

    meta {
        description: "Partition a manifest into chunks."
    }

    parameter_meta {
        manifest: "Manifest to partition."
        N: "Number of chunks to partition into."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        File manifest
        Int N

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        split -a 5 -d --additional-suffix=.txt -e -n l/~{N} ~{manifest} manifest_chunk_
    >>>

    output {
        Array[File] manifest_chunks = glob("manifest_chunk_*.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
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

task GetBasecallModel {
    meta {
        desciption: "Getting the basecall model string of an ONT BAM"
    }
    parameter_meta {
        bam: {
            desciption: "BAM to operate on",
            localization_optional: true
        }
        runid_2_model: "The basecall model for each run."
    }
    input {
        File bam
    }
    output {
        Map[String, String] runid_2_model = read_map("results.tsv")
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep "^@RG" > one_rg_per_line.txt

        while IFS= read -r line
        do
            echo "$line" | tr '\t' '\n' | grep "^DS:" | sed "s/^DS://" | tr ' ' '\n' > tmp.txt
            runid=$(grep "^runid=" tmp.txt | awk -F '=' '{print $2}')
            model=$(grep "^basecall_model=" tmp.txt | awk -F '=' '{print $2}')
            echo -e "${runid}\t${model}" >> results.tsv
        done < one_rg_per_line.txt
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
}
