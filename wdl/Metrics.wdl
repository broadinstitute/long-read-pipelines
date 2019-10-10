version 1.0

import "Structs.wdl"

workflow MetricsOnUnalignedBam {
    input {
        File input_bam
        File input_bai
    }


}

task MakeChrIntervalList {
    input {
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        grep '^@SQ' ~{ref_dict} | awk '{ print $2 "\t" 1 "\t" $3 }' | sed 's/[SL]N://g' | grep -v -e random -e chrUn -e decoy -e alt -e HLA -e EBV > chrs.txt
    >>>

    output {

    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kgarimella/lr-metrics:0.01.00"
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

task ComputeYield {
    input {
        File bam
        File bai

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ~{bam} | awk '{ print $1, length($10) }' | gzip -1 > read_names_and_lengths.txt.gz
    >>>

    output {
        File reads = "reads.fastq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kgarimella/lr-metrics:0.01.00"
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

# yield
# read length histogram (https://www.pacb.com/wp-content/uploads/HiFi-Base-Yield-Dens-plot_Sequel-II-System.png)
# read length N50
# number of zmws
# 3000x3000 zmw yield grid

# depth of coverage (autosome)
# depth of coverage (X)
# depth of coverage (Y)
# depth of coverage (MT)
# read lengths
# alignment stats
# yield
# number of passes
# error rate
# error rate / number of passes
# indel error lengths
# contamination
