version 1.0

import "../../structs/Structs.wdl"

workflow MergeVCFs {

    meta {
        description: "Merges VCFs 'horizonally' (e.g., make a single VCF from many single-sample VCFs)"
    }

    input {
        Array[File] inputArray
        Int chunkSize = 0
        RuntimeAttr? runtime_attr_override_chunks
        RuntimeAttr? runtime_attr_override_merge
        RuntimeAttr? runtime_attr_override_final_merge
    }

    parameter_meta {
        inputArray: "The array of VCFs to merge"
        chunkSize: "Optional number of VCFs to merge in a single task"
        runtime_attr_override_chunks: "Override the default runtime attributes for the chunking task"
        runtime_attr_override_merge: "Override the default runtime attributes for the initial merging tasks"
        runtime_attr_override_final_merge: "Override the default runtime attributes for the final merging task"
    }


    call CreateChunks {
        input:
            inputArray = inputArray,
            chunkSize = chunkSize,
            runtime_attr_override = runtime_attr_override_chunks
    }

    scatter ( inputVCFs in CreateChunks.outputArray ) {
        call MergeVCFsTask {
            input:
                inputVCFs = inputVCFs,
                runtime_attr_override = runtime_attr_override_merge
        }
    }

    call MergeVCFsTask as FinalMerge {
        input:
            inputVCFs = MergeVCFsTask.outputVCF,
            runtime_attr_override = runtime_attr_override_final_merge
    }

    output {
        File outputVCF = FinalMerge.outputVCF
    }
}


task CreateChunks {

    meta {
        description: "Breaks an array into chunks"
    }

    input {
        Array[String] inputArray
        Int chunkSize = 0
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        inputArray: "The array to chunk"
        chunkSize: "Optional chunk size override"
        runtime_attr_override: "Override the default runtime attributes"
    }

    command <<<
        set -euo pipefail
        mkdir chunks
        cd chunks
        lines=`echo "~{length(inputArray)}" | \
                awk 'BEGIN{lines=~{chunkSize}}
                  {if ( !lines ) lines=int(sqrt($1)+.99);
                   while ( $1%lines==1 ) lines += 1;
                   print lines}'`
        split -a4 -l$lines "~{write_lines(inputArray)}"
        for fil in x*;do cat $fil | awk 'BEGIN{sep=""}{printf "%s%s",sep,$0;sep="\t"}END{printf "\n"}' >> ../inputs.tsv; done
    >>>

    output {
        Array[Array[String]] outputArray = read_tsv("inputs.tsv")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "ubuntu:latest"
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

task MergeVCFsTask {

    meta {
        description: "Merges vcfs"
    }

    input {
        Array[File] inputVCFs
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        inputVCFs: {help: "The vcfs to merge -- the vcfs may be given as gs: paths",
                    localization_optional: true}
        runtime_attr_override: "Override the default runtime attributes"
    }

    command <<<
        set -euo pipefail
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools merge --no-index -0 --threads 3 -Oz "~{sep='" "' inputVCFs}" > merged.vcf.gz
    >>>

    output {
        File outputVCF = "merged.vcf.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             4,
        disk_gb:            ceil(size(inputVCFs, "GB")) + 10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/gatk-sv/samtools-cloud:2023-02-01-v0.26.8-beta-9b25c72d"
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
