version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBFindAndSummarizeBAMs {

    meta {
        description: "A workflow that finds BAM files in a directory and emits some basic stats."
    }
    parameter_meta {
        gcs_input_dir:    "GCS path to input BAM files"
        length_threshold: "The minimum length of reads to consider in summary"

        gcs_out_root_dir: "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    input {
        String gcs_input_dir
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBFindAndSummarizeBAMs"

    call Utils.ListFilesOfType {
        input:
            gcs_dir = gcs_input_dir,
            suffixes = [".unaligned.bam"],
            recurse = true
    }

    call SplitFileList {
        input:
            files = ListFilesOfType.files,
            num_sublists = 200
    }

    scatter (filelist in SplitFileList.sublists) {
        call IndexAll { input: filelist = filelist, outdir = outdir }
    }
}

task SplitFileList {
    input {
        Array[String] files
        Int num_sublists

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        FILE_LIST=~{write_lines(files)}

        i=0
        while read file; do
            sublist=$((i % ~{num_sublists}))
            echo $file >> sublist_${sublist}.txt
            i=$((i + 1))
        done < ${FILE_LIST}
    >>>

    output {
        Array[File] sublists = glob("sublist_*.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.40"
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

task IndexAll {
    input {
        File filelist
        String outdir 

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        while read bam; do
            gsutil cp ${bam} .

            bam_basename=$(basename ${bam})
            pbindex ${bam_basename}

            python3 /usr/local/bin/compute_pbi_stats.py -q 0 -l 0 ${bam_basename}.pbi > ${bam_basename}.stats.full.txt
            python3 /usr/local/bin/compute_pbi_stats.py -q 0 -l 10000 ${bam_basename}.pbi > ${bam_basename}.stats.filtered.txt

            gsutil cp ${bam_basename}.stats.full.txt ${outdir}/
            gsutil cp ${bam_basename}.stats.filtered.txt ${outdir}/

            rm ${bam_basename}
        done < ${filelist}
    >>>

    output {
        Array[File] pbis = glob("*.pbi")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            20,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.40"
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
