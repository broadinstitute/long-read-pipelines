version 1.0

##########################################################################################
# A task that runs QUAST to evaluate a given set of assemblies
# on a species with existing reference assembly.
# - Entire Quast output will be tarballed
##########################################################################################

import "../../structs/Structs.wdl"

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

task QuastBenchmark {
    input {
        File? ref
        Array[File] assemblies
        Boolean is_large = false
        Array[String]? labels
        Boolean icarus = false

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        ref:        "reference assembly of the species"
        assemblies: "list of assemblies to evaluate"
        labels:     "Optional list of labels for each assembly (will be used in reports/plots)."
        icarus:     "Include QUAST's Icarus viewer outputs"
    }

    Int minimal_disk_size = 2*(ceil(size(ref, "GB") + size(assemblies, "GB")))
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    String size_optimization = if is_large then "--large" else " "

    command <<<
        set -eux

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        quast ~{true='' false='--no-icarus' icarus} -s \
              "~{size_optimization}" \
              --threads "${num_core}" \
              ~{'-r ' + ref} \
              ~{true='-l ' false='' defined(labels)} ~{sep=', ' labels} \
              ~{sep=' ' assemblies}

        tree -h quast_results/

        if [[ -d quast_results/contigs_reports ]]; then
            tar -zcvf contigs_reports.tar.gz quast_results/contigs_reports
        fi
    >>>

    output {
        File report_txt = "quast_results/latest/report.txt"
        File report_html = "quast_results/latest/report.html"

        Array[File] report_in_various_formats = glob("quast_results/latest/report.*")
        File report_pdf = "quast_results/latest/report.pdf"
        File metrics_tsv = "report_map.txt"

        Array[File] plots = glob("quast_results/latest/basic_stats/*.pdf")
        Array[File] icarus_html = flatten([
            glob("quast_results/latest/icarus.html"),  # glob as workaround for optional output
            glob("quast_results/latest/icarus_viewers/*.html")
        ])
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             16,
        mem_gb:                80,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     0,
        max_retries:           0,
        docker:                "quay.io/biocontainers/quast:5.2.0--py310pl5321hc8f18ef_1"
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

task SummarizeQuastReport {
    input {
        File quast_report_txt
    }

    command <<<
        set -eux
        grep -v -e '^All statistics' -e '^$' ~{quast_report_txt} | \
            sed 's/ /_/g' | \
            sed 's/__\+/\t/g' | \
            sed 's/\s\+$//g' | \
            sed 's/>=/gt/g' | \
            tee report_map.txt

        for i in $(seq 2 $(awk '{print NF}' report_map.txt | sort -nu | tail -n 1))
        do
            j=$(( i - 2 ))  # to make sure the primary, assuming it's the 0-th fed in to this task and the left-most value column
            cut -d$'\t' -f1,${i} < report_map.txt > report_map_${j}.txt
        done
    >>>

    output {
        File quast_metrics_together = "report_map.txt"
        Array[File] quast_metrics = glob("report_map_*.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
