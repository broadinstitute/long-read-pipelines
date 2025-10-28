version 1.0

import "../../structs/Structs.wdl"

task Quast {

    meta {
        description: "A task that runs QUAST to evaluate a given set of assemblies on a species with existing reference assembly. Entire Quast output will be tarballed"
    }
    parameter_meta {
        ref:        "reference assembly of the species"
        assemblies: "list of assemblies to evaluate"
    }

    input {
        File? ref
        Array[File] assemblies
        Boolean is_large = false

        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int minimal_disk_size = 2*(ceil(size(ref, "GB") + size(assemblies, "GB")))
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    String size_optimization = if is_large then "--large" else " "

    command <<<
        set -eux

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        quast --no-icarus \
              "~{size_optimization}" \
              --threads "${num_core}" \
              ~{extra_args} \
              ~{true='-r' false='' defined(ref)} \
              ~{select_first([ref, ""])} \
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

        Array[File] plots = glob("quast_results/latest/basic_stats/*.pdf")

        File? contigs_reports = "contigs_reports.tar.gz"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             16,
        mem_gb:                80,
        disk_gb:               disk_size,
        boot_disk_gb:          25,
        preemptible_tries:     0,
        max_retries:           0,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-quast:5.2.0"
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

    meta {
        description: "A task that summarizes the QUAST report into a single tab-delimited file"
    }
    parameter_meta {
        quast_report_txt: "the QUAST report file"
    }

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
