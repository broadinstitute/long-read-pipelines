version 1.0
import "tasks/Utils.wdl" as Utils

workflow LocalFast5 {
    input {
        Array[String]+ loci
        String aligned_bam
        File   aligned_bai
        String gcs_output_dir
        File summary_txt
    }

    scatter (locus in loci) {
        call Utils.SubsetBam {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                locus = locus
        }
    }

    if (length(loci) > 1)  {
        call Utils.MergeBams {
            input:
                bams = SubsetBam.subset_bam,
                prefix = "merged"
        }
    }

    File subset_bam = select_first([MergeBams.merged_bam, SubsetBam.subset_bam[0]])

    call GetReadnames {
        input:
            bam = subset_bam
    }

    call GetLocalFast5 {
        input:
            readnames = GetReadnames.readnames,
            summary_file = summary_txt,
            gcs_output_dir = gcs_output_dir
    }

}

task GetLocalFast5 {
    input {
        File readnames
        File summary_file
        String gcs_output_dir
    }

    String fast5_dir = sub(summary_file, basename(summary_file), "fast5_pass")

    command <<<
        set -euxo pipefail
        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir fast5
        mkdir output

        while read name; do grep $name ~{summary_file} |cut -f2 >> filenames; done < ~{readnames}
        sort -u filenames -o filenames

        while read filename; do gsutil cp ~{fast5_dir}/$filename fast5/ ; done < filenames

        fast5_subset -i fast5 -s output -l ~{readnames} -t $num_core

        ## save output
        cd output
        gsutil cp *.fast5 ~{gcs_output_dir}
    >>>

    output {
        File local_fast5 = local_fast5
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "nanozoo/ont-fast5-api:3.1.6--a980386"
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

task GetReadnames {
    input {
        File bam
    }

    Int disk_size = 2*size(bam, "GB")

    command <<<
        set -euxo pipefail

        samtools view -f 16 ~{bam} | cut -f1| sort -u > readnames
    >>>

    output {
        File readnames = "readnames"
    }

     #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
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