version 1.0

##########################################################################################
# A task that runs QUAST to evaluate a given set of assemblies
# on a species with existing reference assembly.
# - Entire Quast output will be tarballed
##########################################################################################

import "Structs.wdl"

task Quast {
    input {
        File? ref
        Array[File] assemblies

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        ref:        "reference assembly of the species"
        assemblies: "list of assemblies to evaluate"
    }

    Int disk_size = 2*(ceil(size(ref, "GB") + size(assemblies, "GB")))

    command <<<
        set -x

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        quast --no-icarus \
              --threads $num_core \
              ~{true='-r' false='' defined(ref)} ~{select_first([ref, ""])} ~{sep=' ' assemblies}

        cat quast_results/latest/report.txt | \
            grep -v -e '^All statistics' -e '^$' | \
            sed 's/ /_/g' | \
            sed 's/__\+/\t/g' | \
            sed 's/\s\+$//g' | \
            sed 's/>=/gt/g' | \
            tee report_map.txt
    >>>

    output {
        File report_html = "quast_results/latest/report.html"
        File report_txt = "quast_results/latest/report.txt"
        File report_pdf = "quast_results/latest/report.pdf"
        Array[File] plots = glob("quast_results/latest/basic_stats/*.pdf")

        Map[String, String] metrics = read_map("report_map.txt")
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             2,
        mem_gb:                32,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     0,
        max_retries:           0,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-quast:0.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}
